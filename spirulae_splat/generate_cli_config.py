#!/usr/bin/env python3
"""Generate the standalone CLI trainer's config header from the Python
training config dataclasses.

Source of truth: TrainerConfig (+ presets) in modules/trainer.py and the
nested SpirulaeSplatDataParserConfig / SpirulaeSplatDataManagerConfig /
SpirulaeSplatModelConfig / OptimizerConfig dataclasses. This script parses
them with `ast` (no torch import, so it runs on a fresh checkout before
csrc.so exists) and emits
    spirulae_splat/splat/cuda/csrc/app/generated/cli_config.h
containing:
  - struct SsplatConfig: every training config field, flattened, with the
    Python defaults baked in;
  - SSPLAT_CONFIG_FIELDS(X): an X-macro over (member, cli_key, group,
    choices, help) that the CLI's generic parser/help printer expands;
  - ssplat_apply_preset(): one branch per tyro preset subcommand, assigning
    exactly the fields whose defaults the preset class overrides.

Flattening: nested config names are dropped (--model.sh-degree becomes
--sh-degree). Name collisions across groups must be resolved in RENAMES
below; the script errors on any unlisted collision so a new conflicting
Python field cannot silently shadow an existing flag.

Type mapping:
  int/float/bool/str/Path        -> int / float / bool / std::string
  Optional[int|float]            -> std::optional<int|float>  ("none" on CLI)
  Optional[str|Path]             -> std::string, "" = None
  Literal[str...] (opt. None)    -> std::string + choices ("none" for None)
  Literal[True, False, None]     -> std::optional<bool>
  Tuple[...]                     -> std::array<int|float, N>  (N CLI values)
  Union[bool, int]               -> float via TYPE_OVERRIDES (see entry)
"""

import ast
import os
import sys
from dataclasses import dataclass, field as dc_field
from typing import Optional


REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_PATH = os.path.join(
    REPO, "spirulae_splat/splat/cuda/csrc/app/generated/cli_config.h")

# (group, source file, class name), in flattening order.
CONFIG_SOURCES = [
    ("trainer",     "spirulae_splat/modules/trainer.py",     "TrainerConfig"),
    ("dataparser",  "spirulae_splat/modules/dataparser.py",  "SpirulaeSplatDataParserConfig"),
    ("datamanager", "spirulae_splat/modules/datamanager.py", "SpirulaeSplatDataManagerConfig"),
    ("model",       "spirulae_splat/modules/model.py",       "SpirulaeSplatModelConfig"),
    ("optimizer",   "spirulae_splat/modules/optimizer.py",   "OptimizerConfig"),
]

# TrainerConfig fields that hold nested config objects (flattened separately).
NESTED_CONFIG_FIELDS = {"dataparser", "datamanager", "model", "optimizer"}

# Cross-group name collisions -> flattened CLI/struct name. The script
# errors on collisions not listed here.
RENAMES = {
    # model.split_batch is the engine grad-accumulation path; the
    # datamanager flag is the legacy Python data-path OOM workaround
    # (a no-op with the C++ DataManager).
    ("datamanager", "split_batch"): "dm_split_batch",
}

# Fields whose annotation the generic mapper can't (or shouldn't) handle:
# (group, name) -> (kind, help suffix appended to the docstring).
TYPE_OVERRIDES = {
    # Union[bool, int]: True = probe the image resolution, int = fixed factor.
    ("dataparser", "rescale_camera_to_fit"): (
        ("Float", 0, None),
        " [CLI: 0 = off, -1 = auto-detect from image resolution, > 0 = divide intrinsics by this]"),
}

# Preset subcommands, matching ss_trainer.py's tyro mapping. "3dgs" is the
# base TrainerConfig (no overrides) and the default when no preset is given.
PRESETS = [
    ("3dgs",              "TrainerConfig"),
    ("360-camera",        "TrainerConfig360Camera"),
    ("in-the-wild",       "TrainerConfigInTheWild"),
    ("linear-color",      "TrainerConfigLinear"),
    ("synthetic",         "TrainerConfigSynthetic"),
    ("meshing",           "TrainerConfigMeshing"),
    ("academic-baseline", "TrainerConfigAcademicBaseline"),
]


# ===========================================================================
# AST helpers
# ===========================================================================

class SkipField(Exception):
    pass


REQUIRED = object()


def safe_eval(node):
    """Evaluate a default-value expression without importing the module.
    Supports the literal subset the configs actually use."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, (ast.Tuple, ast.List)):
        return tuple(safe_eval(e) for e in node.elts)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return -safe_eval(node.operand)
    if isinstance(node, ast.BinOp):
        a, b = safe_eval(node.left), safe_eval(node.right)
        if isinstance(node.op, ast.Add):  return a + b
        if isinstance(node.op, ast.Sub):  return a - b
        if isinstance(node.op, ast.Mult): return a * b
        if isinstance(node.op, ast.Div):  return a / b
        if isinstance(node.op, ast.Pow):  return a ** b
    if isinstance(node, ast.Call):
        fname = node.func.id if isinstance(node.func, ast.Name) else None
        if fname == "float" and len(node.args) == 1:
            return float(safe_eval(node.args[0]))
        if fname == "int" and len(node.args) == 1:
            return int(safe_eval(node.args[0]))
        if fname == "Path" and len(node.args) == 1:
            return str(safe_eval(node.args[0]))
        if fname == "field":
            raise SkipField()  # nested config default_factory
    raise ValueError(f"unsupported default expression: {ast.dump(node)}")


@dataclass
class FieldSpec:
    group: str
    pyname: str            # name inside its Python config class
    cname: str = ""        # flattened struct member / CLI key (after renames)
    kind: str = ""         # Int/Float/Bool/String/OptInt/OptFloat/OptBool/TupleI/TupleF
    arity: int = 0         # tuple arity
    choices: tuple = ()    # allowed strings for Literal fields ("" entries excluded)
    allow_none: bool = False   # string field where "none" -> "" is legal
    default: object = None
    required: bool = False
    help: str = ""


def classify_annotation(node, group, name):
    """Return (kind, arity, choices) for an annotation AST node."""
    if (group, name) in TYPE_OVERRIDES:
        return TYPE_OVERRIDES[(group, name)][0]
    if isinstance(node, ast.Name):
        base = {"int": "Int", "float": "Float", "bool": "Bool",
                "str": "String", "Path": "String"}.get(node.id)
        if base:
            return (base, 0, None)
        raise ValueError(f"{group}.{name}: unsupported annotation {node.id}")
    if isinstance(node, ast.Subscript):
        outer = node.value.id if isinstance(node.value, ast.Name) else None
        inner = node.slice
        if outer == "Optional":
            k, _, _ = classify_annotation(inner, group, name)
            if k == "Int":    return ("OptInt", 0, None)
            if k == "Float":  return ("OptFloat", 0, None)
            if k == "String": return ("String", 0, "none")  # "" = None
            raise ValueError(f"{group}.{name}: unsupported Optional[{k}]")
        if outer == "Literal":
            elts = inner.elts if isinstance(inner, ast.Tuple) else [inner]
            vals = [safe_eval(e) for e in elts]
            if set(vals) <= {True, False, None}:
                return ("OptBool", 0, None)
            strs = tuple(v for v in vals if isinstance(v, str))
            has_none = any(v is None for v in vals)
            if len(strs) + has_none == len(vals):
                return ("String", 0, (strs, has_none))
            raise ValueError(f"{group}.{name}: unsupported Literal {vals}")
        if outer == "Tuple":
            elts = inner.elts if isinstance(inner, ast.Tuple) else [inner]
            kinds = [classify_annotation(e, group, name)[0] for e in elts]
            if all(k == "Int" for k in kinds):
                return ("TupleI", len(kinds), None)
            if all(k in ("Int", "Float") for k in kinds):
                return ("TupleF", len(kinds), None)
            raise ValueError(f"{group}.{name}: unsupported Tuple {kinds}")
        if outer == "Union":
            elts = inner.elts if isinstance(inner, ast.Tuple) else [inner]
            names = [e.id if isinstance(e, ast.Name) else None for e in elts]
            raise ValueError(
                f"{group}.{name}: Union[{names}] needs a TYPE_OVERRIDES entry")
    raise ValueError(f"{group}.{name}: unsupported annotation {ast.dump(node)}")


def collapse_doc(doc):
    return " ".join(doc.split())


def find_class(tree, class_name):
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            return node
    raise ValueError(f"class {class_name} not found")


def extract_fields(tree, group, class_name):
    """Ordered FieldSpec list for one config dataclass."""
    cls = find_class(tree, class_name)
    fields = []
    body = cls.body
    for i, stmt in enumerate(body):
        if not (isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name)):
            continue
        name = stmt.target.id
        if group == "trainer" and name in NESTED_CONFIG_FIELDS:
            continue
        # Docstring: string literal expression immediately following.
        doc = ""
        if i + 1 < len(body):
            nxt = body[i + 1]
            if (isinstance(nxt, ast.Expr) and isinstance(nxt.value, ast.Constant)
                    and isinstance(nxt.value.value, str)):
                doc = collapse_doc(nxt.value.value)

        try:
            kind, arity, extra = classify_annotation(stmt.annotation, group, name)
        except SkipField:
            continue
        choices, allow_none = (), False
        if kind == "String" and extra == "none":
            allow_none = True
        elif kind == "String" and isinstance(extra, tuple):
            choices, allow_none = extra

        if (group, name) in TYPE_OVERRIDES:
            doc += TYPE_OVERRIDES[(group, name)][1]

        default, required = None, False
        if stmt.value is None:
            required = True
        else:
            try:
                default = safe_eval(stmt.value)
            except SkipField:
                continue

        fields.append(FieldSpec(
            group=group, pyname=name, kind=kind, arity=arity,
            choices=choices, allow_none=allow_none,
            default=default, required=required, help=doc))
    return fields


def extract_preset_overrides(tree, class_name):
    """{(group, pyname): value} for a preset subclass's default_factory
    keyword overrides. Base TrainerConfig -> {}."""
    if class_name == "TrainerConfig":
        cls = find_class(tree, class_name)
        return {}, collapse_doc(ast.get_docstring(cls) or "")
    cls = find_class(tree, class_name)
    doc = collapse_doc(ast.get_docstring(cls) or "")
    overrides = {}
    for stmt in cls.body:
        if not (isinstance(stmt, ast.AnnAssign) and isinstance(stmt.target, ast.Name)):
            continue
        group = stmt.target.id
        if group not in NESTED_CONFIG_FIELDS:
            raise ValueError(f"{class_name}: unexpected preset field {group}")
        # field(default_factory=lambda: SomeConfig(kw=..., ...))
        call = stmt.value
        assert isinstance(call, ast.Call) and call.func.id == "field", \
            f"{class_name}.{group}: expected field(default_factory=lambda: ...)"
        factory = next(kw.value for kw in call.keywords
                       if kw.arg == "default_factory")
        assert isinstance(factory, ast.Lambda), \
            f"{class_name}.{group}: expected a lambda default_factory"
        inner = factory.body
        assert isinstance(inner, ast.Call), \
            f"{class_name}.{group}: expected a config-constructor call"
        for kw in inner.keywords:
            overrides[(group, kw.arg)] = safe_eval(kw.value)
    return overrides, doc


# ===========================================================================
# C++ emission
# ===========================================================================

def c_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def float_lit(v):
    v = float(v)
    if v != v:
        raise ValueError("NaN default not supported")
    if v == float("inf"):
        return "std::numeric_limits<float>::infinity()"
    if v == float("-inf"):
        return "-std::numeric_limits<float>::infinity()"
    s = repr(v)
    return s + ("f" if ("." in s or "e" in s or "E" in s) else ".0f")


def cpp_value(f: FieldSpec, v):
    if f.kind == "Int":
        return str(int(v))
    if f.kind == "Float":
        # bool default for Union[bool,int] overrides: False -> 0, True -> -1 (auto)
        if isinstance(v, bool):
            v = -1.0 if v else 0.0
        return float_lit(v)
    if f.kind == "Bool":
        return "true" if v else "false"
    if f.kind == "String":
        return c_str("" if v is None else str(v))
    if f.kind == "OptInt":
        return "std::nullopt" if v is None else str(int(v))
    if f.kind == "OptFloat":
        return "std::nullopt" if v is None else float_lit(v)
    if f.kind == "OptBool":
        return "std::nullopt" if v is None else ("true" if v else "false")
    if f.kind == "TupleI":
        return "{" + ", ".join(str(int(x)) for x in v) + "}"
    if f.kind == "TupleF":
        return "{" + ", ".join(float_lit(x) for x in v) + "}"
    raise ValueError(f.kind)


def cpp_type(f: FieldSpec):
    return {
        "Int": "int", "Float": "float", "Bool": "bool",
        "String": "std::string",
        "OptInt": "std::optional<int>",
        "OptFloat": "std::optional<float>",
        "OptBool": "std::optional<bool>",
        "TupleI": f"std::array<int, {f.arity}>",
        "TupleF": f"std::array<float, {f.arity}>",
    }[f.kind]


def choices_str(f: FieldSpec):
    parts = list(f.choices)
    if f.allow_none:
        parts.append("none")
    return "|".join(parts)


def write_if_changed(path, text):
    old = None
    if os.path.exists(path):
        with open(path) as fp:
            old = fp.read()
    if old != text:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fp:
            fp.write(text)
        return True
    return False


def main():
    trees = {}
    all_fields = []
    for group, rel, class_name in CONFIG_SOURCES:
        path = os.path.join(REPO, rel)
        with open(path) as fp:
            trees[group] = ast.parse(fp.read())
        all_fields += extract_fields(trees[group], group, class_name)

    # Flatten with collision detection.
    by_cname = {}
    for f in all_fields:
        f.cname = RENAMES.get((f.group, f.pyname), f.pyname)
        if f.cname in by_cname:
            other = by_cname[f.cname]
            raise SystemExit(
                f"error: flattened name collision '{f.cname}' between "
                f"{other.group}.{other.pyname} and {f.group}.{f.pyname}; "
                f"add a RENAMES entry in {os.path.basename(__file__)}")
        by_cname[f.cname] = f

    # Presets (parsed from trainer.py).
    presets = []
    for preset_name, class_name in PRESETS:
        overrides, doc = extract_preset_overrides(trees["trainer"], class_name)
        items = []
        for (group, pyname), value in overrides.items():
            key = RENAMES.get((group, pyname), pyname)
            if key not in by_cname or by_cname[key].group != group:
                raise SystemExit(
                    f"error: preset {preset_name} overrides unknown field "
                    f"{group}.{pyname}")
            items.append((by_cname[key], value))
        presets.append((preset_name, doc, items))

    # ---- Emit ----
    L = []
    L.append("#pragma once")
    L.append("")
    L.append("// AUTO-GENERATED by spirulae_splat/generate_cli_config.py -- DO NOT EDIT.")
    L.append("// Source of truth: the Python training config dataclasses")
    L.append("// (TrainerConfig + presets, SpirulaeSplatDataParserConfig,")
    L.append("//  SpirulaeSplatDataManagerConfig, SpirulaeSplatModelConfig,")
    L.append("//  OptimizerConfig). Re-run the generator after editing those.")
    L.append("")
    L.append("#include <array>")
    L.append("#include <limits>")
    L.append("#include <optional>")
    L.append("#include <string>")
    L.append("")
    L.append("struct SsplatConfig {")
    cur_group = None
    for f in all_fields:
        if f.group != cur_group:
            cur_group = f.group
            src = dict((g, c) for g, _, c in CONFIG_SOURCES)[cur_group]
            L.append(f"    // ==== {cur_group} ({src}) ====")
        default = "{}" if f.required else cpp_value(f, f.default)
        req = "  // REQUIRED" if f.required else ""
        L.append(f"    {cpp_type(f)} {f.cname} = {default};{req}")
    L.append("};")
    L.append("")
    L.append("// X(member, cli_key, pyname, group, choices, help). cli_key uses '_';")
    L.append("// the parser treats '-' and '_' as equivalent. pyname is the field's")
    L.append("// original name inside its Python config class (differs from cli_key")
    L.append("// only for RENAMES entries); group is the Python sub-config it lives")
    L.append("// in. choices is a '|' list for string fields ('' = free-form);")
    L.append("// 'none' selects the empty string.")
    L.append("#define SSPLAT_CONFIG_FIELDS(X) \\")
    for f in all_fields:
        L.append(f"    X({f.cname}, {c_str(f.cname)}, {c_str(f.pyname)}, "
                 f"{c_str(f.group)}, {c_str(choices_str(f))}, {c_str(f.help)}) \\")
    L.append("    /* end */")
    L.append("")
    L.append("// Required fields (no Python default). Checked after flag parsing.")
    req_names = [f.cname for f in all_fields if f.required]
    L.append("#define SSPLAT_CONFIG_REQUIRED_FIELDS(X) \\")
    for name in req_names:
        L.append(f"    X({name}) \\")
    L.append("    /* end */")
    L.append("")
    L.append("struct SsplatPresetInfo { const char* name; const char* help; };")
    L.append("inline constexpr SsplatPresetInfo kSsplatPresets[] = {")
    for preset_name, doc, _ in presets:
        L.append(f"    {{{c_str(preset_name)}, {c_str(doc)}}},")
    L.append("};")
    L.append("")
    L.append("// Apply a preset's default overrides (tyro subcommand equivalent).")
    L.append("// Returns false for an unknown preset name. \"3dgs\" is the base config.")
    L.append("inline bool ssplat_apply_preset(SsplatConfig& c, const std::string& name) {")
    for preset_name, _, items in presets:
        L.append(f"    if (name == {c_str(preset_name)}) {{")
        for f, value in items:
            L.append(f"        c.{f.cname} = {cpp_value(f, value)};")
        L.append("        return true;")
        L.append("    }")
    L.append("    (void)c;")
    L.append("    return false;")
    L.append("}")
    L.append("")

    changed = write_if_changed(OUT_PATH, "\n".join(L))
    n_presets = len(presets)
    print(f"generate_cli_config: {len(all_fields)} fields, {n_presets} presets"
          f" -> {os.path.relpath(OUT_PATH, REPO)}"
          f" ({'updated' if changed else 'unchanged'})")


if __name__ == "__main__":
    main()
