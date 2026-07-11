// Color-gamut matrices (source gamut -> Rec.709 / sRGB primaries), matching
// spirulae's get_color_transform_matrix (_wrapper_per_pixel.py). Stored
// row-major here; toColMajor() produces the column-major layout GLSL expects
// for `uGamut * rgb`.

export const GAMUTS = {
  'Rec.709':      [1,0,0, 0,1,0, 0,0,1],
  'ACES2065-1':   [ 2.5247180476,-1.1325619434,-0.3921561044,
                   -0.2776344819, 1.3709123773,-0.0932778953,
                   -0.0165202369,-0.1479259606, 1.1644461975],
  'ACEScg':       [ 1.7072552160,-0.6200352595,-0.0872199564,
                   -0.1311566587, 1.1391010566,-0.0079443978,
                   -0.0245499075,-0.1248045805, 1.1493544880],
  'Rec.2020':     [ 1.6604910021,-0.5876411388,-0.0728498633,
                   -0.1245504745, 1.1328998971,-0.0083494226,
                   -0.0181507634,-0.1005788980, 1.1187296614],
  'AdobeRGB':     [ 1.3983671735,-0.3983451225, 0.0000054016,
                   -0.0000103176, 0.9999916496,-0.0000039459,
                   -0.0000003709,-0.0429269510, 1.0429319656],
  'DCI-P3':       [ 1.1548337042,-0.1451763523,-0.0096573518,
                   -0.0393300117, 1.0378282998, 0.0015017119,
                   -0.0184786235,-0.0689101110, 1.0873887345],
};

export function toColMajor(rowMajor) {
  const m = rowMajor;
  return new Float32Array([
    m[0], m[3], m[6],
    m[1], m[4], m[7],
    m[2], m[5], m[8],
  ]);
}

export function hexToRgb(hex) {
  const h = hex.replace('#','');
  return [parseInt(h.slice(0,2),16)/255, parseInt(h.slice(2,4),16)/255, parseInt(h.slice(4,6),16)/255];
}
