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


#define USE_ANYHIT 0

int BoundingBoxIntersectCount = 0;
vec2 _BoundingBoxIntersect(
	vec3 minCorner, vec3 maxCorner, vec3 rayOrigin, vec3 invDir
) {
	BoundingBoxIntersectCount += 1;
	vec3 near = (minCorner - rayOrigin) * invDir;
	vec3 far  = (maxCorner - rayOrigin) * invDir;
	
	vec3 tmin = min(near, far);
	vec3 tmax = max(near, far);
	
	float t0 = max( max(tmin.x, tmin.y), tmin.z);
	float t1 = min( min(tmax.x, tmax.y), tmax.z);
	
	return max(t0, 0.0) > t1 ? vec2(INFINITY) : vec2(t0, t1);
	// return max(t0, 0.0) > t1 ? INFINITY : t0;
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



uvec4 GetBoxNodeData(const in float i) {
	float ix2 = i;
	ivec2 uv0 = ivec2( mod(ix2 + 0.0, AABB_TEXTURE_SIZE), (ix2 + 0.0) / AABB_TEXTURE_SIZE );
	return texelFetch(tAABBTexture, uv0, 0);
}

vec2 BoundingBoxIntersect(
	uvec4 nodeData, vec3 rayOrigin, vec3 invDir
) {
	vec2 x = unpackHalf2x16(nodeData.y), y = unpackHalf2x16(nodeData.z), z = unpackHalf2x16(nodeData.w);
	vec3 pmin = vec3(x.x, y.x, z.x);
	vec3 pmax = vec3(x.y, y.y, z.y);
	return _BoundingBoxIntersect(pmin, pmax, rayOrigin, invDir);
}



#if USE_ANYHIT

struct StackData {
	float ptr;
	float t0, t1;
};

StackData createStackData(float ptr, uvec4 currentBoxNodeData, vec3 rayOrigin, vec3 inverseDir) {
	vec2 t = BoundingBoxIntersect(currentBoxNodeData, rayOrigin, inverseDir);
	return StackData(ptr, t.x, t.y);
}

#else  // USE_ANYHIT

struct StackData {
	float ptr;
	float t0;
};

StackData createStackData(float ptr, uvec4 currentBoxNodeData, vec3 rayOrigin, vec3 inverseDir) {
	return StackData(ptr, BoundingBoxIntersect(currentBoxNodeData, rayOrigin, inverseDir).x);
}

#endif  // USE_ANYHIT


// Ray-scene intersection

float stackLevels[32];

float SceneIntersect() {
	uvec4 currentBoxNodeData, nodeAData, nodeBData, tmpNodeData;
	
	vec3 inverseDir = 1.0 / rayDirection;

	StackData currentStackData, stackDataA, stackDataB, tmpStackData;

	float d, t = INFINITY;
	float stackptr = 0.0;
	float splatID = -1.0;
	vec4 splatRGBA = vec4(0);

	int skip = FALSE;
	int splatLookupNeeded = FALSE;

	currentBoxNodeData = GetBoxNodeData(stackptr);
	currentStackData = createStackData(stackptr, currentBoxNodeData, rayOrigin, inverseDir);
	stackLevels[0] = currentStackData.ptr;
	skip = (currentStackData.t0 < t) ? TRUE : FALSE;

	while (true) {
		if (skip == FALSE) {
			// decrease pointer by 1 (0.0 is root level)
			if (--stackptr < 0.0) // went past the root level, terminate loop
				break;

			currentStackData.ptr = stackLevels[int(stackptr)];

			currentBoxNodeData = GetBoxNodeData(currentStackData.ptr);
			currentStackData = createStackData(currentStackData.ptr, currentBoxNodeData, rayOrigin, inverseDir);

			if (currentStackData.t0 >= t)
				continue;
		}
		skip = FALSE; // reset skip

		// negative is child, positive is triangle
		float id = float(int(currentBoxNodeData.x));

		if (id < 0.0) // < 0.0 signifies an inner node
		{
			nodeAData = GetBoxNodeData(currentStackData.ptr + 1.0);
			nodeBData = GetBoxNodeData(-id);
			stackDataA = createStackData(currentStackData.ptr + 1.0, nodeAData, rayOrigin, inverseDir);
			stackDataB = createStackData(-id, nodeBData, rayOrigin, inverseDir);

			// first sort the branch node data so that A is the closest
			if (stackDataB.t0 < stackDataA.t0)
			{
				tmpStackData = stackDataB, stackDataB = stackDataA, stackDataA = tmpStackData;
				tmpNodeData = nodeBData, nodeBData = nodeAData, nodeAData = tmpNodeData;
			} // branch 'b' now has the larger rayT value of A and B

			if (stackDataB.t0 < t) // see if branch B (the larger rayT) needs to be processed
			{
				currentStackData = stackDataB;
				currentBoxNodeData = nodeBData;
				skip = TRUE; // this will prevent the stackptr from decreasing by 1
			}
			if (stackDataA.t0 < t) // see if branch A (the smaller rayT) needs to be processed 
			{
				if (skip == TRUE) // if larger branch B needed to be processed also,
					stackLevels[int(stackptr++)] = stackDataB.ptr; // cue larger branch B for future round
				
				currentStackData = stackDataA;
				currentBoxNodeData = nodeAData;
				skip = TRUE; // this will prevent the stackptr from decreasing by 1
			}

			continue;
		} // end if (id < 0.0) // inner node

		// else this is a leaf

		if (id != splatID) {

			ivec2 uv0 = ivec2( mod(2.0*id+0.0, TRI_TEXTURE_SIZE), (2.0*id+0.0) / TRI_TEXTURE_SIZE );
			ivec2 uv1 = ivec2( mod(2.0*id+1.0, TRI_TEXTURE_SIZE), (2.0*id+1.0) / TRI_TEXTURE_SIZE );
			uvec4 vd0 = texelFetch(tTriangleTexture, uv0, 0);
			uvec4 vd1 = texelFetch(tTriangleTexture, uv1, 0);

			vec3 pos = vec3(uintBitsToFloat(vd0.x), uintBitsToFloat(vd0.y), uintBitsToFloat(vd0.z));
			vec3 au = vec3(unpackHalf2x16(vd0.w), unpackHalf2x16(vd1.x).x);
			vec3 av = vec3(unpackHalf2x16(vd1.x).y, unpackHalf2x16(vd1.y));
			vec2 uv;
			d = SplatIntersect(pos, au, av, rayOrigin, rayDirection, uv);

			if (d < t) {
				t = d;
				splatID = id;
				splatRGBA = vec4(unpackHalf2x16(vd1.z), unpackHalf2x16(vd1.w));
				splatRGBA.w *= max(1.0-dot(uv, uv), 0.0);
				splatLookupNeeded = TRUE;
			}
		}

	} // end while (TRUE)

	if (splatLookupNeeded == TRUE)
	{
		vec4 rgba = splatRGBA;
		hitColor = rgba.xyz;
		hitOpacity = rgba.w;
	}

	return t;
}


// Ray-scene intersection, any hit, assuming ID don't repeat

#if USE_ANYHIT

#define N_INT 32
vec2 intersects[N_INT];  // (id, depth)
vec4 intersectRGBA[N_INT];
int intersectPtr = 0;

void SceneIntersectAnyHit() {
	uvec4 currentBoxNodeData, nodeAData, nodeBData, tmpNodeData;
	
	vec3 inverseDir = 1.0 / rayDirection;

	StackData currentStackData, stackDataA, stackDataB, tmpStackData;

	float d, t = INFINITY;
	float stackptr = 0.0;
	intersectPtr = 0;

	int skip = FALSE;

	currentBoxNodeData = GetBoxNodeData(stackptr);
	currentStackData = createStackData(stackptr, currentBoxNodeData, rayOrigin, inverseDir);
	stackLevels[0] = currentStackData.ptr;
	skip = (currentStackData.t0 < t) ? TRUE : FALSE;

	while (true) {
		if (skip == FALSE) {
			// decrease pointer by 1 (0.0 is root level)
			if (--stackptr < 0.0) // went past the root level, terminate loop
				break;

			currentStackData.ptr = stackLevels[int(stackptr)];

			currentBoxNodeData = GetBoxNodeData(currentStackData.ptr);
			currentStackData = createStackData(currentStackData.ptr, currentBoxNodeData, rayOrigin, inverseDir);

			if (currentStackData.t0 >= t)
				continue;
		}
		skip = FALSE; // reset skip

		// negative is child, positive is triangle
		float id = float(int(currentBoxNodeData.x));

		if (id < 0.0) // < 0.0 signifies an inner node
		{
			nodeAData = GetBoxNodeData(currentStackData.ptr + 1.0);
			nodeBData = GetBoxNodeData(-id);
			stackDataA = createStackData(currentStackData.ptr + 1.0, nodeAData, rayOrigin, inverseDir);
			stackDataB = createStackData(-id, nodeBData, rayOrigin, inverseDir);

			// first sort the branch node data so that A is the closest
			// if (stackDataB.t0 < stackDataA.t0)
			if (stackDataB.t1 < stackDataA.t1)
			{
				tmpStackData = stackDataB, stackDataB = stackDataA, stackDataA = tmpStackData;
				tmpNodeData = nodeBData, nodeBData = nodeAData, nodeAData = tmpNodeData;
			} // branch 'b' now has the larger rayT value of A and B

			if (stackDataB.t0 < t) // see if branch B (the larger rayT) needs to be processed
			{
				currentStackData = stackDataB;
				currentBoxNodeData = nodeBData;
				skip = TRUE; // this will prevent the stackptr from decreasing by 1
			}
			if (stackDataA.t0 < t) // see if branch A (the smaller rayT) needs to be processed 
			{
				if (skip == TRUE) // if larger branch B needed to be processed also,
					stackLevels[int(stackptr++)] = stackDataB.ptr; // cue larger branch B for future round
				
				currentStackData = stackDataA;
				currentBoxNodeData = nodeAData;
				skip = TRUE; // this will prevent the stackptr from decreasing by 1
			}

			continue;
		} // end if (id < 0.0) // inner node

		// else this is a leaf

		ivec2 uv0 = ivec2( mod(2.0*id+0.0, TRI_TEXTURE_SIZE), (2.0*id+0.0) / TRI_TEXTURE_SIZE );
		ivec2 uv1 = ivec2( mod(2.0*id+1.0, TRI_TEXTURE_SIZE), (2.0*id+1.0) / TRI_TEXTURE_SIZE );
		uvec4 vd0 = texelFetch(tTriangleTexture, uv0, 0);
		uvec4 vd1 = texelFetch(tTriangleTexture, uv1, 0);

		vec3 pos = vec3(uintBitsToFloat(vd0.x), uintBitsToFloat(vd0.y), uintBitsToFloat(vd0.z));
		vec3 au = vec3(unpackHalf2x16(vd0.w), unpackHalf2x16(vd1.x).x);
		vec3 av = vec3(unpackHalf2x16(vd1.x).y, unpackHalf2x16(vd1.y));
		vec2 uv;
		d = SplatIntersect(pos, au, av, rayOrigin, rayDirection, uv);

		if (d < t) {
			// t = d;
			vec4 rgba = vec4(unpackHalf2x16(vd1.z), unpackHalf2x16(vd1.w));
			rgba.w *= max(1.0-dot(uv, uv), 0.0);
			if (rgba.w > 0.0) {
				// insertion sort
				vec2 iz = vec2(intersectPtr, d);
				int j = intersectPtr-1;
				#if 0  // do it later seems to be faster
				while (j >= 0 && intersects[j].y > d) {
					intersects[j+1] = intersects[j];
					j--;
				}
				#endif
				intersects[j+1] = iz;

				intersectRGBA[intersectPtr] = rgba;
				if (++intersectPtr >= N_INT)
					return;
			}
		}

	} // end while (TRUE)

}

#endif


vec4 CalculateRadiance() {
	vec3 cumColor = vec3(0);
	float T = 1.0;

	// prevent unrolling, save compile time
	int ZERO = max(-int(focal.x), 0);

#if !USE_ANYHIT
	// ray trace

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

#else
	// any hit + sort

	SceneIntersectAnyHit();

	for (int i = ZERO+1; i < intersectPtr; ++i) {
        vec2 iz = intersects[i];
        int j = i-1;

        while (j >= 0 && intersects[j].y > iz.y) {
            intersects[j+1] = intersects[j];
            j--;
        }
        intersects[j+1] = iz;
    }

	for (int i = ZERO; i < intersectPtr; i++) {
		vec4 rgba = intersectRGBA[int(intersects[i].x)];
		cumColor += T * rgba.w * rgba.xyz;
		T *= (1.0-rgba.w);
	}

#endif

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
