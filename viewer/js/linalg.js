// Minimal, dependency-free linear algebra for the viewer.
// Vectors are plain arrays; matrices are column-major Float32Array-friendly
// arrays (mat4 = 16, mat3 = 9), matching WebGL uniform layout.

export const v3 = {
  add: (a, b) => [a[0]+b[0], a[1]+b[1], a[2]+b[2]],
  sub: (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]],
  scale: (a, s) => [a[0]*s, a[1]*s, a[2]*s],
  dot: (a, b) => a[0]*b[0]+a[1]*b[1]+a[2]*b[2],
  cross: (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]],
  len: (a) => Math.hypot(a[0], a[1], a[2]),
  norm: (a) => { const l = Math.hypot(a[0],a[1],a[2]) || 1; return [a[0]/l, a[1]/l, a[2]/l]; },
  lerp: (a, b, t) => [a[0]+(b[0]-a[0])*t, a[1]+(b[1]-a[1])*t, a[2]+(b[2]-a[2])*t],
};

// Quaternions stored as [x, y, z, w].
export const quat = {
  identity: () => [0, 0, 0, 1],
  norm: (q) => { const l = Math.hypot(q[0],q[1],q[2],q[3]) || 1; return [q[0]/l,q[1]/l,q[2]/l,q[3]/l]; },
  mul: (a, b) => [
    a[3]*b[0] + a[0]*b[3] + a[1]*b[2] - a[2]*b[1],
    a[3]*b[1] - a[0]*b[2] + a[1]*b[3] + a[2]*b[0],
    a[3]*b[2] + a[0]*b[1] - a[1]*b[0] + a[2]*b[3],
    a[3]*b[3] - a[0]*b[0] - a[1]*b[1] - a[2]*b[2],
  ],
  fromAxisAngle: (axis, angle) => {
    const a = v3.norm(axis), s = Math.sin(angle/2);
    return [a[0]*s, a[1]*s, a[2]*s, Math.cos(angle/2)];
  },
  rotVec: (q, v) => {
    // v + 2*q.w*(qv x v) + 2*(qv x (qv x v))
    const qv = [q[0], q[1], q[2]];
    const t = v3.scale(v3.cross(qv, v), 2);
    return v3.add(v3.add(v, v3.scale(t, q[3])), v3.cross(qv, t));
  },
  // 3x3 rotation matrix (column-major) from quaternion (cam-to-world)
  toMat3: (q) => {
    const [x,y,z,w] = q;
    const xx=x*x, yy=y*y, zz=z*z, xy=x*y, xz=x*z, yz=y*z, wx=w*x, wy=w*y, wz=w*z;
    return [
      1-2*(yy+zz), 2*(xy+wz),   2*(xz-wy),
      2*(xy-wz),   1-2*(xx+zz), 2*(yz+wx),
      2*(xz+wy),   2*(yz-wx),   1-2*(xx+yy),
    ];
  },
  fromMat3: (m) => {
    // m column-major
    const m00=m[0], m10=m[1], m20=m[2], m01=m[3], m11=m[4], m21=m[5], m02=m[6], m12=m[7], m22=m[8];
    const tr = m00+m11+m22;
    let x,y,z,w;
    if (tr > 0) { const s = Math.sqrt(tr+1)*2; w=0.25*s; x=(m21-m12)/s; y=(m02-m20)/s; z=(m10-m01)/s; }
    else if (m00>m11 && m00>m22) { const s=Math.sqrt(1+m00-m11-m22)*2; w=(m21-m12)/s; x=0.25*s; y=(m01+m10)/s; z=(m02+m20)/s; }
    else if (m11>m22) { const s=Math.sqrt(1+m11-m00-m22)*2; w=(m02-m20)/s; x=(m01+m10)/s; y=0.25*s; z=(m12+m21)/s; }
    else { const s=Math.sqrt(1+m22-m00-m11)*2; w=(m10-m01)/s; x=(m02+m20)/s; y=(m12+m21)/s; z=0.25*s; }
    return quat.norm([x,y,z,w]);
  },
};

export const mat3 = {
  transpose: (m) => [m[0],m[3],m[6], m[1],m[4],m[7], m[2],m[5],m[8]],
  mul: (a, b) => {
    const o = new Array(9);
    for (let c=0;c<3;c++) for (let r=0;r<3;r++)
      o[c*3+r] = a[r]*b[c*3] + a[3+r]*b[c*3+1] + a[6+r]*b[c*3+2];
    return o;
  },
  mulVec: (m, v) => [
    m[0]*v[0]+m[3]*v[1]+m[6]*v[2],
    m[1]*v[0]+m[4]*v[1]+m[7]*v[2],
    m[2]*v[0]+m[5]*v[1]+m[8]*v[2],
  ],
};

export const mat4 = {
  identity: () => new Float32Array([1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1]),
  mul: (a, b) => {
    const o = new Float32Array(16);
    for (let c=0;c<4;c++) for (let r=0;r<4;r++) {
      let s=0; for (let k=0;k<4;k++) s += a[k*4+r]*b[c*4+k];
      o[c*4+r]=s;
    }
    return o;
  },
  // world->cam from rotation (cam-to-world, mat3 col-major) and camera position
  view: (R, pos) => {
    // view = [R^T | -R^T pos]
    const Rt = mat3.transpose(R);
    const t = mat3.mulVec(Rt, pos);
    return new Float32Array([
      Rt[0], Rt[1], Rt[2], 0,
      Rt[3], Rt[4], Rt[5], 0,
      Rt[6], Rt[7], Rt[8], 0,
      -t[0], -t[1], -t[2], 1,
    ]);
  },
  perspective: (fovy, aspect, near, far) => {
    const f = 1/Math.tan(fovy/2);
    const nf = 1/(near-far);
    return new Float32Array([
      f/aspect,0,0,0,
      0,f,0,0,
      0,0,(far+near)*nf,-1,
      0,0,2*far*near*nf,0,
    ]);
  },
  ortho: (l,r,b,t,near,far) => {
    const lr=1/(l-r), bt=1/(b-t), nf=1/(near-far);
    return new Float32Array([
      -2*lr,0,0,0,
      0,-2*bt,0,0,
      0,0,2*nf,0,
      (l+r)*lr,(t+b)*bt,(far+near)*nf,1,
    ]);
  },
};
