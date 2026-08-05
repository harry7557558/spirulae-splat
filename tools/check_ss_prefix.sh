#!/bin/bash
# Fail if this project defines an SS_* name a system header already owns.
#
# SS_ is a short prefix: <signal.h> and winuser.h both have names in it, and a
# collision surfaces as a redefinition warning at best and as silently wrong
# code at worst. The reserved set lives in tools/ss_reserved_names.txt.
#
# Only *declarations* are checked -- `#define SS_X`, `option(SS_X ...)`,
# `set(SS_X ...)`. Using a Windows style constant is fine; defining one is not.
#
# Usage:  bash tools/check_ss_prefix.sh

cd "$(dirname "$0")/.." || exit 1

reserved=$(grep -oE '^SS_[A-Z0-9_]+' tools/ss_reserved_names.txt | sort -u)
[ -z "$reserved" ] && { echo "no reserved names loaded" >&2; exit 1; }

declared=$(git grep -hIoE \
    -e '^[[:space:]]*#[[:space:]]*define[[:space:]]+SS_[A-Za-z0-9_]+' \
    -e '(option|set)\([[:space:]]*SS_[A-Za-z0-9_]+' \
    -- ':!tools/ss_reserved_names.txt' ':!tools/check_ss_prefix.sh' \
    | grep -oE 'SS_[A-Za-z0-9_]+' | sort -u)

hits=$(comm -12 <(echo "$reserved") <(echo "$declared"))

if [ -n "$hits" ]; then
    echo "Reserved SS_* names defined by this project:"
    echo ""
    echo "$hits" | sed 's/^/  /'
    echo ""
    echo "These belong to <signal.h> or winuser.h. Pick another name; see"
    echo "tools/ss_reserved_names.txt."
    exit 1
fi

echo "OK: no reserved SS_* names defined ($(echo "$declared" | wc -l) checked)."
