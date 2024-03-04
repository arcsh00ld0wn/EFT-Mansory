#pragma once
// Libraries
#include <math.h>
#include <assert.h>

class FVector;
class FRotator;
class Vector3;

float DegToRad(float x);
float RadToDeg(float x);
float DistancePointToLine(FVector Point, FVector LineOrigin, FVector Dir);

float Vec3Dot(const Vector3* pV1, const Vector3* pV2);

// Class Vector2
// This class represents a 2D vector.
class Vector2 {

public:

	// -------------------- Attributes -------------------- //

	// Components of the vector
	float x, y;


	// -------------------- Methods -------------------- //

	// Constructor
	Vector2(float x = 0, float y = 0) : x(x), y(y) {}

	// Constructor
	Vector2(const Vector2& vector) : x(vector.x), y(vector.y) {}

	// + operator
	Vector2 operator+(const Vector2& v) const {
		return Vector2(x + v.x, y + v.y);
	}

	// += operator
	Vector2& operator+=(const Vector2& v) {
		x += v.x; y += v.y;
		return *this;
	}

	// - operator
	Vector2 operator-(const Vector2& v) const {
		return Vector2(x - v.x, y - v.y);
	}

	// -= operator
	Vector2& operator-=(const Vector2& v) {
		x -= v.x; y -= v.y;
		return *this;
	}

	// = operator
	Vector2& operator=(const Vector2& vector) {
		if (&vector != this) {
			x = vector.x;
			y = vector.y;
		}
		return *this;
	}

	// == operator
	bool operator==(const Vector2& v) const {
		return x == v.x && y == v.y;
	}

	// * operator
	Vector2 operator*(float f) const {
		return Vector2(f * x, f * y);
	}

	// *= operator
	Vector2& operator*=(float f) {
		x *= f; y *= f;
		return *this;
	}

	// / operator
	Vector2 operator/(float f) const {
		assert(f != 0);
		float inv = 1.f / f;
		return Vector2(x * inv, y * inv);
	}

	// /= operator
	Vector2& operator/=(float f) {
		assert(f != 0);
		float inv = 1.f / f;
		x *= inv; y *= inv;
		return *this;
	}

	// - operator
	Vector2 operator-() const {
		return Vector2(-x, -y);
	}

	// [] operator
	float& operator[](int i) {
		assert(i >= 0 && i <= 1);
		switch (i) {
		case 0: return x;
		case 1: return y;
		}
		return y;
	}

	// Normalize the vector and return it
	Vector2 normalize() {
		float l = length();
		assert(l > 0);
		x /= l;
		y /= l;
		return *this;
	}

	// Clamp the vector values between 0 and 1
	Vector2 clamp01() {
		if (x > 1.f) x = 1.f;
		else if (x < 0.f) x = 0.f;
		if (y > 1.f) y = 1.f;
		else if (y < 0.f) y = 0.f;
		return *this;
	}

	// Return the squared length of the vector
	float lengthSquared() const { return x * x + y * y; }

	// Return the length of the vector
	float length() const { return sqrt(lengthSquared()); }
};


//Vector3
class Vector3
{
public:
	Vector3() : x(0.f), y(0.f), z(0.f)
	{

	}

	Vector3(float _x, float _y, float _z) : x(_x), y(_y), z(_z)
	{

	}
	~Vector3()
	{

	}

	float x;
	float y;
	float z;

	float Dot(Vector3 v);
	Vector3 operator+(Vector3 v);
	Vector3 operator*(const float other);
	Vector3 operator-(Vector3 v);
	float Distance(Vector3 v);

	float Area()
	{
		return x * x + y * y + z * z;
	}

	float Magnitude()
	{
		return sqrt(Area());
	}

	Vector3& Normalize()
	{
		float norm = Magnitude();
		Vector3 c = { x / norm, y / norm, z / norm };
		return c;
	}

	Vector3& Multiply(Vector3 v)
	{
		Vector3 c = { x * v.x, y * v.y, z * v.z };
		return c;
	}

	Vector3& Multiply(float f)
	{
		Vector3 c = { x * f, y * f, z * f };
		return c;
	}

	Vector3& Divide(float d)
	{
		Vector3 c = { x / d, y / d, z / d };
		return c;
	}

	Vector3& Subtract(Vector3 source)
	{
		Vector3 c = { x - source.x, y - source.y, z - source.z };
		return c;
	}



};

class Vector4
{
public:
	float x, z, y, w;

	Vector4();
	Vector4(float x, float y, float z, float w);
	Vector4(float x, float y, float z);
	Vector4(float* x, float* y, float* z, float* w);
	~Vector4() {}

	Vector4 operator= (const Vector4& v);
	const Vector4 operator* (const float& scalar) const;
	const Vector4 operator+ (const Vector4& v) const;

	//inline void show() { std::cout << this->x << " " << this->y << " " << this->z << " " << this->w << std::endl; }

private:

	// only to allow do that: m(0) = Vector4(1,1,1,1) //
	float* px;
	float* py;
	float* pz;
	float* pw;

	bool pointer;  // to check what constructor was called //
};

class FVector
{
public:
	float x, y, z;

	FVector();
	FVector(float x, float y, float z);
	FVector(const FVector& other);
	FVector(const Vector3 vec);

	FVector operator+ (const FVector& other) const;
	FVector operator- (const FVector& other) const;
	FVector operator* (const float other) const;
	float operator* (const FVector& other) const;

	bool operator == (const FVector& other) const;
	bool operator != (const FVector& other) const;

	FVector& operator= (const FVector& other);
	FVector& operator+= (const FVector& other);
	FVector& operator-= (const FVector& other);
	FVector& operator*= (const float other);

	float& operator[](size_t i);
	const float& operator[](size_t i) const;

	float GetLength() const;
	float Dot(FVector v);
	float Distance(FVector v);
	float Distance2(FVector a, FVector b);
	float GetMagnitudeSqr();
	FVector normalize();

	FRotator VectorAngles() const;
};

class FMatrix
{
public:
	FMatrix() : m{
		{ 0.f, 0.f, 0.f, 0.f },
		{ 0.f, 0.f, 0.f, 0.f },
		{ 0.f, 0.f, 0.f, 0.f },
		{ 0.f, 0.f, 0.f, 0.f } }
	{
	}

	FMatrix(const FMatrix&) = default;


	float* operator[](size_t i) { return m[i]; }
	const float* operator[](size_t i) const { return m[i]; }


	FVector operator*(const FVector& vec);
	FMatrix operator*(const FMatrix& other);
	float m[4][4];
};

class FRotator
{
public:
	float pitch;
	float yaw;
	float roll;

	FRotator();
	FRotator(float pitch, float yaw, float roll);
	FRotator(const FRotator& other);

	void ToSourceAngles();
	void ToUnityAngles();
	void Normalize();
	FVector AngleVector();
	void AngleVectors(FVector* x, FVector* y, FVector* z);
};

struct FQuat
{
	float x;
	float y;
	float z;
	float w;

	FQuat operator*(const FQuat& other);
};

struct FTransform
{
public:
	FQuat Rotation;
	FVector Translation;
private:
	float pad0;
public:
	FVector Scale3D;
private:
	float pad1;

public:
	FMatrix ToMatrixWithScale();

};

struct FBoxSphereBounds
{
	FVector	Origin;
	FVector BoxExtent;
	float SphereRadius;
};

#pragma once

#include <sstream>

class Vector
{
public:
	Vector(void)
	{
		Invalidate();
	}
	Vector(float X, float Y, float Z)
	{
		x = X;
		y = Y;
		z = Z;
	}
	Vector(const float* clr)
	{
		x = clr[0];
		y = clr[1];
		z = clr[2];
	}

	void Init(float ix = 0.0f, float iy = 0.0f, float iz = 0.0f)
	{
		x = ix; y = iy; z = iz;
	}
	bool IsValid() const
	{
		return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
	}
	void Invalidate()
	{
		x = y = z = std::numeric_limits<float>::infinity();
	}

	float& operator[](int i)
	{
		return ((float*)this)[i];
	}
	float operator[](int i) const
	{
		return ((float*)this)[i];
	}

	void Zero()
	{
		x = y = z = 0.0f;
	}

	bool operator==(const Vector& src) const
	{
		return (src.x == x) && (src.y == y) && (src.z == z);
	}
	bool operator!=(const Vector& src) const
	{
		return (src.x != x) || (src.y != y) || (src.z != z);
	}
	inline float Distance(const Vector& vector)
	{
		return sqrt(
			(x - vector.x) * (x - vector.x) +
			(y - vector.y) * (y - vector.y) +
			(z - vector.z) * (z - vector.z));
	}
	Vector& operator+=(const Vector& v)
	{
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Vector& operator-=(const Vector& v)
	{
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	Vector& operator*=(float fl)
	{
		x *= fl;
		y *= fl;
		z *= fl;
		return *this;
	}
	Vector& operator*=(const Vector& v)
	{
		x *= v.x;
		y *= v.y;
		z *= v.z;
		return *this;
	}
	Vector& operator/=(const Vector& v)
	{
		x /= v.x;
		y /= v.y;
		z /= v.z;
		return *this;
	}
	Vector& operator+=(float fl)
	{
		x += fl;
		y += fl;
		z += fl;
		return *this;
	}
	Vector& operator/=(float fl)
	{
		x /= fl;
		y /= fl;
		z /= fl;
		return *this;
	}
	Vector& operator-=(float fl)
	{
		x -= fl;
		y -= fl;
		z -= fl;
		return *this;
	}

	void Clamp()
	{
		if (this->x > 180.f)
			this->x -= 360.f;

		else if (this->x < -180.f)
			this->x += 360.f;

		if (this->z > 180.f)
			this->z -= 360.f;

		else if (this->z < -180.f)
			this->z += 360.f;

		if (this->x < -89.f)
			this->x = -89.f;

		if (this->x > 89.f)
			this->x = 89.f;

		while (this->z < -180.0f)
			this->z += 360.0f;

		while (this->z > 180.0f)
			this->z -= 360.0f;
	}

	void NormalizeInPlace()
	{
		*this = Normalized();
	}
	Vector Normalized() const
	{
		Vector res = *this;
		float l = res.Length();
		if (l != 0.0f) {
			res /= l;
		}
		else {
			res.x = res.y = res.z = 0.0f;
		}
		return res;
	}

	float DistTo(const Vector& vOther) const
	{
		Vector delta;

		delta.x = x - vOther.x;
		delta.y = y - vOther.y;
		delta.z = z - vOther.z;

		return delta.Length();
	}
	float DistToSqr(const Vector& vOther) const
	{
		Vector delta;

		delta.x = x - vOther.x;
		delta.y = y - vOther.y;
		delta.z = z - vOther.z;

		return delta.LengthSqr();
	}
	float Dot(const Vector& vOther) const
	{
		return (x * vOther.x + y * vOther.y + z * vOther.z);
	}
	float Length() const
	{
		return sqrt(x * x + y * y + z * z);
	}
	float LengthSqr(void) const
	{
		return (x * x + y * y + z * z);
	}
	float Length2D() const
	{
		return sqrt(x * x + y * y);
	}

	Vector& operator=(const Vector& vOther)
	{
		x = vOther.x; y = vOther.y; z = vOther.z;
		return *this;
	}

	Vector operator-(void) const
	{
		return Vector(-x, -y, -z);
	}
	Vector operator+(const Vector& v) const
	{
		return Vector(x + v.x, y + v.y, z + v.z);
	}
	Vector operator-(const Vector& v) const
	{
		return Vector(x - v.x, y - v.y, z - v.z);
	}
	Vector operator*(float fl) const
	{
		return Vector(x * fl, y * fl, z * fl);
	}
	Vector operator*(const Vector& v) const
	{
		return Vector(x * v.x, y * v.y, z * v.z);
	}
	Vector operator/(float fl) const
	{
		return Vector(x / fl, y / fl, z / fl);
	}
	Vector operator/(const Vector& v) const
	{
		return Vector(x / v.x, y / v.y, z / v.z);
	}
	inline Vector Normalize()
	{
		Vector vector;
		float length = this->Length();

		if (length != 0) {
			vector.x = x / length;
			vector.y = y / length;
			vector.z = z / length;
		}
		else
			vector.x = vector.y = 0.0f; vector.z = 1.0f;

		return vector;
	}
	inline float Normalizes()
	{
		Vector res = *this;
		float l = res.Length();
		if (l != 0.0f)
		{
			res /= l;
		}
		else
		{
			res.x = res.y = res.z = 0.0f;
		}
		return l;
	}
	float x, y, z;
};

inline Vector operator*(float lhs, const Vector& rhs)
{
	return rhs * lhs;
}
inline Vector operator/(float lhs, const Vector& rhs)
{
	return rhs / lhs;
}

class __declspec(align(16)) VectorAligned : public Vector
{
public:
	inline VectorAligned(void) {};
	inline VectorAligned(float X, float Y, float Z)
	{
		Init(X, Y, Z);
	}

public:
	explicit VectorAligned(const Vector& vOther)
	{
		Init(vOther.x, vOther.y, vOther.z);
	}

	VectorAligned& operator=(const Vector& vOther)
	{
		Init(vOther.x, vOther.y, vOther.z);
		return *this;
	}

	VectorAligned& operator=(const VectorAligned& vOther)
	{
		Init(vOther.x, vOther.y, vOther.z);
		return *this;
	}

	float w;
};

class vector2
{
public:
	vector2()
	{
		x = y = 0.f;
	}

	vector2(float fx, float fy)
	{
		x = fx;
		y = fy;
	}

	float x, y;

	vector2 operator+(const vector2& input) const
	{
		return vector2{ x + input.x, y + input.y };
	}

	vector2 operator-(const vector2& input) const
	{
		return vector2{ x - input.x, y - input.y };
	}

	vector2 operator/(float input) const
	{
		return vector2{ x / input, y / input };
	}

	vector2 operator*(float input) const
	{
		return vector2{ x * input, y * input };
	}

	vector2& operator-=(const vector2& v)
	{
		x -= v.x;
		y -= v.y;
		return *this;
	}

	vector2& operator/=(float input)
	{
		x /= input;
		y /= input;
		return *this;
	}

	vector2& operator*=(float input)
	{
		x *= input;
		y *= input;
		return *this;
	}

	float length() const
	{
		return std::sqrt((x * x) + (y * y));
	}

	vector2 normalized() const
	{
		return { x / length(), y / length() };
	}

	float dot_product(vector2 input) const
	{
		return (x * input.x) + (y * input.y);
	}

	float distance(vector2 input) const
	{
		return (*this - input).length();
	}

	bool empty() const
	{
		return x == 0.f && y == 0.f;
	}
};

class vector3
{
public:
	vector3()
	{
		x = y = z = 0.f;
	}

	vector3(float fx, float fy, float fz)
	{
		x = fx;
		y = fy;
		z = fz;
	}

	float x, y, z;

	vector3 operator+(const vector3& input) const
	{
		return vector3{ x + input.x, y + input.y, z + input.z };
	}

	vector3 operator-(const vector3& input) const
	{
		return vector3{ x - input.x, y - input.y, z - input.z };
	}

	vector3 operator/(float input) const
	{
		return vector3{ x / input, y / input, z / input };
	}

	vector3 operator*(float input) const
	{
		return vector3{ x * input, y * input, z * input };
	}

	vector3& operator-=(const vector3& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;

		return *this;
	}

	vector3& operator/=(float input)
	{
		x /= input;
		y /= input;
		z /= input;
		return *this;
	}

	vector3& operator*=(float input)
	{
		x *= input;
		y *= input;
		z *= input;
		return *this;
	}

	bool operator==(const vector3& input) const
	{
		return x == input.x && y == input.y && z == input.z;
	}

	void make_absolute()
	{
		x = std::abs(x);
		y = std::abs(y);
		z = std::abs(z);
	}

	float length_sqr() const
	{
		return (x * x) + (y * y) + (z * z);
	}

	float length() const
	{
		return std::sqrt(length_sqr());
	}

	float length_2d() const
	{
		return std::sqrt((x * x) + (y * y));
	}

	vector3 normalized() const
	{
		return { x / length(), y / length(), z / length() };
	}

	float dot_product(vector3 input) const
	{
		return (x * input.x) + (y * input.y) + (z * input.z);
	}

	float distance(vector3 input) const
	{
		return (*this - input).length();
	}

	float distance_2d(vector3 input) const
	{
		return (*this - input).length_2d();
	}

	void clamp()
	{
		while (x > 89.0f)
			x -= 180.0f;

		while (x < -89.0f)
			x += 180.0f;

		while (y > 180.0f)
			y -= 360.0f;

		while (y < -180.0f)
			y += 360.0f;

		z = 0.f;
	}

	bool empty() const
	{
		return x == 0.f && y == 0.f && z == 0.f;
	}
};

class vector4
{
public:
	float x;
	float y;
	float z;
	float w;

	vector4();
	vector4(float x, float y, float z, float w);

	vector4 operator+(const vector4& vector) const;
	vector4 operator-(const vector4& vector) const;
	vector4 operator-() const;
	vector4 operator*(float number) const;
	vector4 operator/(float number) const;

	vector4& operator+=(const vector4& vector);
	vector4& operator-=(const vector4& vector);
	vector4& operator*=(float number);
	vector4& operator/=(float number);

	bool operator==(const vector4& vector) const;
	bool operator!=(const vector4& vector) const;

	inline float Dot(const vector4& vector)
	{
		return x * vector.x + y * vector.y + z * vector.z + w * vector.w;
	}

	inline float Distance(const vector4& vector)
	{
		return sqrt(
			(x - vector.x) * (x - vector.x) +
			(y - vector.y) * (y - vector.y) +
			(z - vector.z) * (z - vector.z) +
			(w - vector.w) * (w - vector.w));
	}

	bool empty() const
	{
		return x == 0.f && y == 0.f && z == 0.f && w == 0.f;
	}
};

inline bool vector4::operator==(const vector4& vector) const
{
	return x == vector.x && y == vector.y && z == vector.z && w == vector.w;
}

inline bool vector4::operator!=(const vector4& vector) const
{
	return x != vector.x || y != vector.y || z != vector.z || w != vector.w;
}

inline vector4 vector4::operator+(const vector4& vector) const
{
	return vector4(x + vector.x, y + vector.y, z + vector.z, w + vector.w);
}

inline vector4 vector4::operator-(const vector4& vector) const
{
	return vector4(x - vector.x, y - vector.y, z - vector.z, w - vector.w);
}

inline vector4 vector4::operator-() const
{
	return vector4(-x, -y, -z, -w);
}

inline vector4 vector4::operator*(float number) const
{
	return vector4(x * number, y * number, z * number, w * number);
}

inline vector4 vector4::operator/(float number) const
{
	return vector4(x / number, y / number, z / number, w / number);
}

class matrix
{
public:
	inline float* operator[](int i)
	{
		return m[i];
	}

	inline const float* operator[](int i) const
	{
		return m[i];
	}

	inline float* Base()
	{
		return &m[0][0];
	}

	inline const float* Base() const
	{
		return &m[0][0];
	}
public:

	inline matrix()
	{
	}

	inline matrix(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33)
	{
		Init(
			m00, m01, m02, m03,
			m10, m11, m12, m13,
			m20, m21, m22, m23,
			m30, m31, m32, m33
		);
	}

	inline void Init(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
	)
	{
		m[0][0] = m00;
		m[0][1] = m01;
		m[0][2] = m02;
		m[0][3] = m03;

		m[1][0] = m10;
		m[1][1] = m11;
		m[1][2] = m12;
		m[1][3] = m13;

		m[2][0] = m20;
		m[2][1] = m21;
		m[2][2] = m22;
		m[2][3] = m23;

		m[3][0] = m30;
		m[3][1] = m31;
		m[3][2] = m32;
		m[3][3] = m33;
	}

	matrix transpose() const
	{
		return matrix(
			m[0][0], m[1][0], m[2][0], m[3][0],
			m[0][1], m[1][1], m[2][1], m[3][1],
			m[0][2], m[1][2], m[2][2], m[3][2],
			m[0][3], m[1][3], m[2][3], m[3][3]);
	}

	float m[4][4];
};

using vec2_t = vector2;
using vec3_t = vector3;
using vec4_t = vector4;
using mat4x4_t = matrix;

struct TransformAccessReadOnly
{
	uintptr_t pTransformData;
};

struct TransformData
{
	uintptr_t pTransformArray;
	uintptr_t pTransformIndices;
};

struct Matrix34
{
	vector4 vec0;
	vector4 vec1;
	vector4 vec2;
};