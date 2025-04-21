#version 300 es

precision highp float;
precision highp int;

#include "shader-utils.glsl"
#line 8

uniform vec2 focal;
uniform vec2 viewport;
uniform int camera_model;
uniform vec4 distortion;
uniform vec4 undistortion;
uniform mat4 projection, view;

out vec4 fragColor;

//-----------------------------------------------------------------------

// Based on https://github.com/erichlof/THREE.js-PathTracing-Renderer

#define INFINITY 1048576.0
#define TRUE 1
#define FALSE 0

uniform highp usampler2D tTriangleTexture;
uniform highp usampler2D tAABBTexture;

#define TRI_TEXTURE_SIZE 2048.0
#define AABB_TEXTURE_SIZE 4096.0

vec3 rayOrigin, rayDirection;
vec3 hitColor;
float hitOpacity;


int BoundingBoxIntersectCount = 0;
float BoundingBoxIntersect(
	vec3 minCorner, vec3 maxCorner, vec3 rayOrigin, vec3 invDir
) {
	BoundingBoxIntersectCount += 1;
	vec3 near = (minCorner - rayOrigin) * invDir;
	vec3 far  = (maxCorner - rayOrigin) * invDir;
	
	vec3 tmin = min(near, far);
	vec3 tmax = max(near, far);
	
	float t0 = max( max(tmin.x, tmin.y), tmin.z);
	float t1 = min( min(tmax.x, tmax.y), tmax.z);
	
	return max(t0, 0.0) > t1 ? INFINITY : t0;
	//return max(t0, 0.0) <= t1 ? t0 : INFINITY;
}


int BVH_TriangleIntersectCount = 0;
float BVH_TriangleIntersect(
	vec3 v0, vec3 v1, vec3 v2,
	vec3 rayOrigin, vec3 rayDirection,
	out float u, out float v
) {
	BVH_TriangleIntersectCount += 1;
	vec3 edge1 = v1 - v0;
	vec3 edge2 = v2 - v0;
	vec3 pvec = cross(rayDirection, edge2);
	float det = 1.0 / dot(edge1, pvec);
	vec3 tvec = rayOrigin - v0;
	u = dot(tvec, pvec) * det;
	vec3 qvec = cross(tvec, edge1);
	v = dot(rayDirection, qvec) * det;
	float t = dot(edge2, qvec) * det;
	return (t <= 0.0 || u < 0.0 || u > 1.0 || v < 0.0 || u + v > 1.0) ? INFINITY : t;
}


int SplatIntersectCount = 0;
float SplatIntersect(
	vec3 pos, vec3 au, vec3 av,
	vec3 ro, vec3 rd,
    out vec2 uv
) {
	SplatIntersectCount += 1;
	// ro + rd * t = pos + au * u + av * v
	// [au av rd] * [u v -t]' = ro - pos
    vec3 uvt = inverse(mat3(au, av, rd)) * (ro - pos);
    uv = uvt.xy;
    float t = -uvt.z;
	return (t <= 0.0 || dot(uv, uv) >= 1.0) ? INFINITY : t;
}


#if 0
// uncompressed tAABBTexture sampler
void GetBoxNodeData(const in float i, inout vec4 boxNodeData0, inout vec4 boxNodeData1) {
	// each bounding box's data is encoded in 2 rgba(or xyzw) texture slots
	float ix2 = i * 2.0;
	// (ix2 + 0.0) corresponds to .x = idTriangle,  .y = aabbMin.x, .z = aabbMin.y, .w = aabbMin.z 
	// (ix2 + 1.0) corresponds to .x = idRightChild .y = aabbMax.x, .z = aabbMax.y, .w = aabbMax.z 

	ivec2 uv0 = ivec2( mod(ix2 + 0.0, AABB_TEXTURE_SIZE), (ix2 + 0.0) / AABB_TEXTURE_SIZE ); // data0
	ivec2 uv1 = ivec2( mod(ix2 + 1.0, AABB_TEXTURE_SIZE), (ix2 + 1.0) / AABB_TEXTURE_SIZE ); // data1

	boxNodeData0 = texelFetch(tAABBTexture, uv0, 0);
	boxNodeData1 = texelFetch(tAABBTexture, uv1, 0);
}
#else
// compressed tAABBTexture usampler, only 28fps -> 29fps ??
void GetBoxNodeData(const in float i, inout vec4 boxNodeData0, inout vec4 boxNodeData1) {
	float ix2 = i;
	ivec2 uv0 = ivec2( mod(ix2 + 0.0, AABB_TEXTURE_SIZE), (ix2 + 0.0) / AABB_TEXTURE_SIZE );
	uvec4 data = texelFetch(tAABBTexture, uv0, 0);

	int idx = int(data.x);
	boxNodeData0.x = float(idx);
	boxNodeData1.x = float(-idx);

	vec2 x = unpackHalf2x16(data.y), y = unpackHalf2x16(data.z), z = unpackHalf2x16(data.w);
	boxNodeData0.yzw = vec3(x.x, y.x, z.x);
	boxNodeData1.yzw = vec3(x.y, y.y, z.y);
}
#endif


// TODO: try stackless one
// Also try do everything along one ray instead of doing multiple bounces

vec2 stackLevels[24];

float SceneIntersect() {
	vec4 currentBoxNodeData0, nodeAData0, nodeBData0, tmpNodeData0;
	vec4 currentBoxNodeData1, nodeAData1, nodeBData1, tmpNodeData1;
	
	vec3 inverseDir = 1.0 / rayDirection;

	vec2 currentStackData, stackDataA, stackDataB, tmpStackData;

	float d, t = INFINITY;
	float stackptr = 0.0;
	float id = 0.0;
	float triangleID = -1.0;
	vec4 triangleRGBA = vec4(0);
	
	int skip = FALSE;
	int triangleLookupNeeded = FALSE;

	GetBoxNodeData(stackptr, currentBoxNodeData0, currentBoxNodeData1);
	currentStackData = vec2(stackptr, BoundingBoxIntersect(currentBoxNodeData0.yzw, currentBoxNodeData1.yzw, rayOrigin, inverseDir));
	stackLevels[0] = currentStackData;
	skip = (currentStackData.y < t) ? TRUE : FALSE;

	while (true) {
		if (skip == FALSE) {
			// decrease pointer by 1 (0.0 is root level)
			if (--stackptr < 0.0) // went past the root level, terminate loop
				break;

			currentStackData = stackLevels[int(stackptr)];
			
			if (currentStackData.y >= t)
				continue;
			
			GetBoxNodeData(currentStackData.x, currentBoxNodeData0, currentBoxNodeData1);
		}
		skip = FALSE; // reset skip

		if (currentBoxNodeData0.x < 0.0) // < 0.0 signifies an inner node
		{
			GetBoxNodeData(currentStackData.x + 1.0, nodeAData0, nodeAData1);
			GetBoxNodeData(currentBoxNodeData1.x, nodeBData0, nodeBData1);
			stackDataA = vec2(currentStackData.x + 1.0, BoundingBoxIntersect(nodeAData0.yzw, nodeAData1.yzw, rayOrigin, inverseDir));
			stackDataB = vec2(currentBoxNodeData1.x, BoundingBoxIntersect(nodeBData0.yzw, nodeBData1.yzw, rayOrigin, inverseDir));

			// first sort the branch node data so that 'a' is the smallest
			if (stackDataB.y < stackDataA.y)
			{
				tmpStackData = stackDataB;
				stackDataB = stackDataA;
				stackDataA = tmpStackData;

				tmpNodeData0 = nodeBData0;   tmpNodeData1 = nodeBData1;
				nodeBData0   = nodeAData0;   nodeBData1   = nodeAData1;
				nodeAData0   = tmpNodeData0; nodeAData1   = tmpNodeData1;
			} // branch 'b' now has the larger rayT value of 'a' and 'b'

			if (stackDataB.y < t) // see if branch 'b' (the larger rayT) needs to be processed
			{
				currentStackData = stackDataB;
				currentBoxNodeData0 = nodeBData0;
				currentBoxNodeData1 = nodeBData1;
				skip = TRUE; // this will prevent the stackptr from decreasing by 1
			}
			if (stackDataA.y < t) // see if branch 'a' (the smaller rayT) needs to be processed 
			{
				if (skip == TRUE) // if larger branch 'b' needed to be processed also,
					stackLevels[int(stackptr++)] = stackDataB; // cue larger branch 'b' for future round
							// also, increase pointer by 1
				
				currentStackData = stackDataA;
				currentBoxNodeData0 = nodeAData0;
				currentBoxNodeData1 = nodeAData1;
				skip = TRUE; // this will prevent the stackptr from decreasing by 1
			}

			continue;
		} // end if (currentBoxNodeData0.x < 0.0) // inner node

		// else this is a leaf

		id = 2.0 * currentBoxNodeData0.x;
		if (id != triangleID) {

			ivec2 uv0 = ivec2( mod(id + 0.0, TRI_TEXTURE_SIZE), (id + 0.0) / TRI_TEXTURE_SIZE );
			ivec2 uv1 = ivec2( mod(id + 1.0, TRI_TEXTURE_SIZE), (id + 1.0) / TRI_TEXTURE_SIZE );
			uvec4 vd0 = texelFetch(tTriangleTexture, uv0, 0);
			uvec4 vd1 = texelFetch(tTriangleTexture, uv1, 0);

			vec3 pos = vec3(uintBitsToFloat(vd0.x), uintBitsToFloat(vd0.y), uintBitsToFloat(vd0.z));
			vec3 au = vec3(unpackHalf2x16(vd0.w), unpackHalf2x16(vd1.x).x);
			vec3 av = vec3(unpackHalf2x16(vd1.x).y, unpackHalf2x16(vd1.y));
			vec2 uv;
			d = SplatIntersect(pos, au, av, rayOrigin, rayDirection, uv);

			if (d < t) {
				t = d;
				triangleID = id;
				triangleRGBA = vec4(unpackHalf2x16(vd1.z), unpackHalf2x16(vd1.w));
				triangleRGBA.w *= max(1.0-dot(uv, uv), 0.0);
				triangleLookupNeeded = TRUE;
			}
		}

	} // end while (TRUE)

	if (triangleLookupNeeded == TRUE)
	{
		vec4 rgba = triangleRGBA;
		hitColor = rgba.xyz;
		hitOpacity = rgba.w;
	}

	return t;
}


vec4 CalculateRadiance() {
	vec3 cumColor = vec3(0);
	float T = 1.0;

	// prevent unrolling, save compile time
	int ZERO = max(-int(focal.x), 0);
	for (int bounces = ZERO; bounces < 8; bounces++)  // 16
	{
		float t = SceneIntersect();
		if (t == INFINITY)
			break;

		cumColor += T * hitOpacity * hitColor;
		T *= 1.0-hitOpacity;
		rayOrigin += (t+1e-6) * rayDirection;

		if (T < 0.05)  // 0.01
			break;
	}

	return vec4(cumColor, 1.0-T);
}



void main( void )
{
	vec3 cameraPosition = vec3(inverse(view)[3]);

    vec2 pos2d = vec2(1,-1)*(gl_FragCoord.xy-0.5*viewport)/focal;
	vec2 pos2d_undist = pos2d;

    float fr = -1.0;
    if (camera_model == 0) {
        fr = opencv_radius(distortion);
        pos2d_undist = undistort_opencv_iterative(pos2d, distortion);
    }
    else if (camera_model == 1) {
        fr = fisheye_radius(0.5*PI, distortion);
        pos2d_undist = undistort_fisheye(pos2d, undistortion);
    }
    const vec3 vignetting_background = vec3(0.1);
    if (fr > 0.0 && length(pos2d) > fr) {
        fragColor = vec4(vignetting_background, 1.0);
        return;
    }

	rayOrigin = cameraPosition;
	rayDirection = transpose(mat3(view)) * normalize(vec3(pos2d_undist, 1));

	vec4 color = CalculateRadiance();
	// color = vec4(vec3(BoundingBoxIntersectCount) / 4096.0, 1.0);
	// color = vec4(pow(vec3(BoundingBoxIntersectCount-1) / 4096.0, vec3(0.3)), 1.0);
	// color.xyz = vec3(SplatIntersectCount) / 1024.0;

    if (fr > 0.0) {
        float vignetting = min(length(pos2d)/fr, 1.0);
        vignetting = pow(1.0-pow(vignetting,20.0), 2.0);
        color.xyz = mix(vignetting_background*color.w, color.xyz, vignetting);
    }

	fragColor = color;
}
