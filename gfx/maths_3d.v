module gfx

import math
import rand


pub struct Point {
pub:
    x f64 @[required]
    y f64 @[required]
    z f64 @[required]
}

pub struct Vector {
pub:
    x f64 @[required]
    y f64 @[required]
    z f64 @[required]
}

// NOTE: take care when directly creating a new Direction
//       the following code assumes x^2+y^2+z^2=1
pub struct Direction {
pub:
    x f64 @[required]
    y f64 @[required]
    z f64 @[required]
}

// NOTE: take care when directly creating a new Normal
//       the following code assumes x^2+y^2+z^2=1
pub struct Normal {
pub:
    x f64 @[required]
    y f64 @[required]
    z f64 @[required]
}

pub struct Ray {
pub:
    e Point      @[required]
    d Direction  @[required]
    t_min f64 = 10e-5        // some small epsilon
    t_max f64 = math.inf(1)  // positive inf
    time f64 = 0.0 
}

pub struct Frame {
pub:
    o Point     = point_origin
    x Direction = direction_x
    y Direction = direction_y
    z Direction = direction_z

    // alternative "LookAt" definition
    // unused if eye == target or up is bad direction (zero length)
    eye    Point     = point_origin
    target Point     = point_origin
    up     Direction = direction_no

    velocity Vector = Vector{ 0.0, 0.0, 0.0 }
}



///////////////////////////////////////////////////////////
// constructors and constants

pub fn direction_xyz(x f64, y f64, z f64) Direction {
    l := math.sqrt(x*x + y*y + z*z)
    return Direction{ x / l, y / l, z / l }
}

pub fn direction_uniform_sphere() Direction {
    mut u := f64(0)
    mut v := f64(0)
    mut l2 := f64(2)
    for l2 > 1 {
        u = 2 * rand.f64() - 1
        v = 2 * rand.f64() - 1
        l2 = u * u + v * v
    }
    s := math.sqrt(1 - l2)
    return Direction{
        x: 2 * u * s,
        y: 2 * v * s,
        z: 1 - 2 * l2,
    }
}

pub fn normal_xyz(x f64, y f64, z f64) Normal {
    l := math.sqrt(x*x + y*y + z*z)
    return Normal{ x / l, y / l, z / l }
}

pub fn frame_oz(o Point, z Direction) Frame {
    v := direction_uniform_sphere()
    y := v.cross(z).direction()
    x := y.cross(z).direction()
    return Frame{ o:o, x:x, y:y, z:z }
}

pub fn frame_oxy(o Point, x Direction, y Direction) Frame {
    z := x.cross(y).direction()
    return Frame{ o:o, x:x, y:y, z:z }
}

pub fn frame_lookat(eye Point, target Point, up Direction) Frame {
    o := eye
    z := target.direction_to(eye)  // point away from scene
    x := up.cross(z).direction()
    y := z.cross(x).direction()
    return Frame{ o:o, x:x, y:y, z:z }
}



pub const (
    point_origin = Point{ 0, 0, 0 }

    direction_x  = Direction{ 1, 0, 0 }
    direction_y  = Direction{ 0, 1, 0 }
    direction_z  = Direction{ 0, 0, 1 }
    direction_no = Direction{ 0, 0, 0 }

    frame_standard = Frame{
        o: point_origin,
        x: direction_x,
        y: direction_y,
        z: direction_z,
    }
)


///////////////////////////////////////////////////////////
// printing methods

pub fn (v Point) str() string {
    return 'Point{ $v.x, $v.y, $v.z }'
}

pub fn (v Vector) str() string {
    return 'Vector{ $v.x, $v.y, $v.z }'
}

pub fn (d Direction) str() string {
    return 'Direction{ $d.x, $d.y, $d.z }'
}

pub fn (d Normal) str() string {
    return 'Normal{ $d.x, $d.y, $d.z }'
}


///////////////////////////////////////////////////////////
// direction and magnitude/length methods

pub fn (v Vector) direction() Direction {
    return direction_xyz(v.x, v.y, v.z)
}

pub fn (v Vector) magnitude() f64 {
    return v.length()
}
pub fn (d Direction) magnitude() f64 {
    return 1
}
pub fn (n Normal) magnitude() f64 {
    return 1
}

pub fn (v Vector) length_squared() f64 {
    // same as: return v.dot(v)
    return v.x * v.x + v.y * v.y + v.z * v.z
}
pub fn (d Direction) length_squared() f64 {
    return 1
}
pub fn (n Normal) length_squared() f64 {
    return 1
}

pub fn (v Vector) length() f64 {
    return math.sqrt(v.length_squared())
}
pub fn (d Direction) length() f64 {
    return 1
}
pub fn (n Normal) length() f64 {
    return 1
}

pub fn (v Vector) l1_norm() f64 {
    return math.abs(v.x) + math.abs(v.y) + math.abs(v.z)
}
pub fn (d Direction) l1_norm() f64 {
    return math.abs(d.x) + math.abs(d.y) + math.abs(d.z)
}
pub fn (n Normal) l1_norm() f64 {
    return math.abs(n.x) + math.abs(n.y) + math.abs(n.z)
}

pub fn (v Vector) l2_norm() f64 {
    return v.length()
}
pub fn (d Direction) l2_norm() f64 {
    return 1
}
pub fn (n Normal) l2_norm() f64 {
    return 1
}

pub fn (v Vector) linf_norm() f64 {
    return max3(math.abs(v.x), math.abs(v.y), math.abs(v.z))
}
pub fn (d Direction) linf_norm() f64 {
    return max3(math.abs(d.x), math.abs(d.y), math.abs(d.z))
}
pub fn (n Normal) linf_norm() f64 {
    return max3(math.abs(n.x), math.abs(n.y), math.abs(n.z))
}

pub fn (v Vector) normalize() Vector {
    l := v.length()
    return Vector{ v.x / l, v.y / l, v.z / l }
}


///////////////////////////////////////////////////////////
// conversion methods

pub fn (v Vector) as_direction() Direction {
    return v.direction()
}
pub fn (v Vector) as_normal() Normal {
    return normal_xyz(v.x, v.y, v.z)
}
pub fn (d Direction) as_vector() Vector {
    return Vector{ d.x, d.y, d.z }
}
pub fn (d Direction) as_normal() Normal {
    return Normal{ d.x, d.y, d.z }
}
pub fn (n Normal) as_direction() Direction {
    return Direction{ n.x, n.y, n.z }
}
pub fn (n Normal) as_vector() Vector {
    return Vector{ n.x, n.y, n.z }
}

pub fn (f Frame) as_frame() Frame {
    if f.eye != f.target && f.up != direction_no {
        return frame_lookat(f.eye, f.target, f.up)
    }
    return f
}


///////////////////////////////////////////////////////////
// convenience methods

pub fn (a Point) vector_to(b Point) Vector {
    return Vector{ b.x - a.x, b.y - a.y, b.z - a.z }
}
pub fn (a Point) direction_to(b Point) Direction {
    return a.vector_to(b).direction()
}
pub fn (a Point) distance_to(b Point) f64 {
    return math.sqrt( a.distance_squared_to(b) )
}
pub fn (a Point) distance_squared_to(b Point) f64 {
    return (
        (a.x - b.x) * (a.x - b.x) +
        (a.y - b.y) * (a.y - b.y) +
        (a.z - b.z) * (a.z - b.z)
    )
}

pub fn (a Point) ray_to(b Point) Ray {
    return Ray{
        e: a,
        d: a.direction_to(b),
        t_max: a.distance_to(b),
    }
}

pub fn (a Point) ray_through(b Point) Ray {
    return Ray{ e: a, d: a.direction_to(b) }
}

pub fn (a Point) ray_along(d Direction) Ray {
    return Ray{ e: a, d: d }
}




///////////////////////////////////////////////////////////
// interpolation methods

pub fn average_point(ps ...Point) Point {
    mut x := f64(0)
    mut y := f64(0)
    mut z := f64(0)
    for p in ps {
        x += p.x
        y += p.y
        z += p.z
    }
    return Point{ x / ps.len, y / ps.len, z / ps.len }
}
pub fn (a Point) average(b Point) Point {
    return Point{
        (a.x + b.x) / 2,
        (a.y + b.y) / 2,
        (a.z + b.z) / 2,
    }
}

pub fn (a Vector) average(b Vector) Vector {
    return Vector{ (a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2 }
}

pub fn (a Direction) average(b Direction) Direction {
    return direction_xyz(a.x + b.x, a.y + b.y, a.z + b.z)
}
pub fn (a Normal) average(b Normal) Normal {
    return normal_xyz(a.x + b.x, a.y + b.y, a.z + b.z)
}

pub fn (from Point) lerp(to Point, factor f64) Point {
    return Point{
        from.x * (1.0 - factor) + to.x * factor,
        from.y * (1.0 - factor) + to.y * factor,
        from.z * (1.0 - factor) + to.z * factor,
    }
}

pub fn (from Vector) lerp(to Vector, factor f64) Vector {
    return Vector{
        from.x * (1.0 - factor) + to.x * factor,
        from.y * (1.0 - factor) + to.y * factor,
        from.z * (1.0 - factor) + to.z * factor,
    }
}


///////////////////////////////////////////////////////////
// arithmetic methods

pub fn (p Point) add(v Vector) Point {
    return Point{ p.x + v.x, p.y + v.y, p.z + v.z }
}
pub fn (a Vector) add(b Vector) Vector {
    return Vector{ a.x + b.x, a.y + b.y, a.z + b.z }
}
pub fn (a Vector) + (b Vector) Vector {
    return Vector{ a.x + b.x, a.y + b.y, a.z + b.z }
}

pub fn (p Point) sub(v Vector) Point {
    return Point{ p.x - v.x, p.y - v.y, p.z - v.z }
}
pub fn (a Vector) sub(b Vector) Vector {
    return Vector{ a.x - b.x, a.y - b.y, a.z - b.z }
}
pub fn (a Vector) - (b Vector) Vector {
    return Vector{ a.x - b.x, a.y - b.y, a.z - b.z }
}

pub fn (a Vector) dot<T>(b T) f64 {
    return a.x * b.x + a.y * b.y + a.z * b.z
}
pub fn (a Direction) dot<T>(b T) f64 {
    return a.x * b.x + a.y * b.y + a.z * b.z
}
pub fn (a Normal) dot<T>(b T) f64 {
    return a.x * b.x + a.y * b.y + a.z * b.z
}

pub fn (a Vector) cross<T>(b T) Vector {
    return Vector{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    }
}
pub fn (a Direction) cross<T>(b T) Vector {
    return Vector{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    }
}
pub fn (a Normal) cross<T>(b T) Vector {
    return Vector{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    }
}


pub fn (a Vector) mult<T>(b T) Vector {
    return Vector{ a.x * b.x, a.y * b.y, a.z * b.z }
}

pub fn (a Vector) div<T>(b T) Vector {
    return Vector{ a.x / b.x, a.y / b.y, a.z / b.z }
}

pub fn (a Vector) scale(s f64) Vector {
    return Vector{ a.x * s, a.y * s, a.z * s }
}
pub fn (d Direction) scale(s f64) Vector {
    return Vector{ d.x * s, d.y * s, d.z * s }
}
pub fn (n Normal) scale(s f64) Vector {
    return Vector{ n.x * s, n.y * s, n.z * s }
}

pub fn (v Vector) negate() Vector {
    return Vector{ -v.x, -v.y, -v.z }
}
pub fn (d Direction) negate() Direction {
    return Direction{ -d.x, -d.y, -d.z }
}
pub fn (n Normal) negate() Normal {
    return Normal{ -n.x, -n.y, -n.z }
}



///////////////////////////////////////////////////////////
// evaluation and t-related methods

pub fn (r Ray) at(t f64) Point {
    return r.e.add(r.d.scale(r.clamp(t)))
}

pub fn (r Ray) clamp(t f64) f64 {
    return math.max(r.t_min, math.min(r.t_max, t))
}

pub fn (r Ray) valid_t(t f64) bool {
    return r.t_min <= t && t <= r.t_max
}


///////////////////////////////////////////////////////////
// transformation methods

// transform point from world to local space
pub fn (f Frame) w2l_point(p Point) Point {
    v := f.o.vector_to(p)
    x := f.x.as_vector().dot(v)
    y := f.y.as_vector().dot(v)
    z := f.z.as_vector().dot(v)
    return Point{ x, y, z }
}
// transform point from local to world space
pub fn (f Frame) l2w_point(p Point) Point {
    x := f.x.scale(p.x)
    y := f.y.scale(p.y)
    z := f.z.scale(p.z)
    return f.o.add(x + y + z)
}

// transform vector from world to local space
pub fn (f Frame) w2l_vector(v Vector) Vector {
    x := f.x.as_vector().dot(v)
    y := f.y.as_vector().dot(v)
    z := f.z.as_vector().dot(v)
    return Vector{ x, y, z }
}
// transform vector from local to world space
pub fn (f Frame) l2w_vector(v Vector) Vector {
    x := f.x.scale(v.x)
    y := f.y.scale(v.y)
    z := f.z.scale(v.z)
    return x + y + z
}

// NOTE: since frame axes are orthonormal and the
//       direction/normal input is normalized, we
//       do not need to normalize the result.

// transform direction from world to local space
pub fn (f Frame) w2l_direction(d Direction) Direction {
    v := d.as_vector()
    x := f.x.as_vector().dot(v)
    y := f.y.as_vector().dot(v)
    z := f.z.as_vector().dot(v)
    return Direction{ x, y, z }
}
// transform direction from local to world space
pub fn (f Frame) l2w_direction(d Direction) Direction {
    vl := d.as_vector()
    vw := f.x.scale(vl.x) + f.y.scale(vl.y) + f.z.scale(vl.z)
    return Direction{ vw.x, vw.y, vw.z }
}

// transform normal from world to local space
pub fn (f Frame) w2l_normal(n Normal) Normal {
    v := n.as_vector()
    x := f.x.as_vector().dot(v)
    y := f.y.as_vector().dot(v)
    z := f.z.as_vector().dot(v)
    return Normal{ x, y, z }
}
// transform normal from local to world space
pub fn (f Frame) l2w_normal(n Normal) Normal {
    vl := n.as_vector()
    vw := f.x.scale(vl.x) + f.y.scale(vl.y) + f.z.scale(vl.z)
    return Normal{ vw.x, vw.y, vw.z }
}

// NOTE: since frame axes are orthonormal and the
//       ray's direction is normalized, we do not
//       need to scale t_min or t_max.

// transform ray from world to local space
pub fn (f Frame) w2l_ray(r Ray) Ray {
    return Ray{
        e: f.w2l_point(r.e),
        d: f.w2l_direction(r.d),
        t_min: r.t_min,
        t_max: r.t_max,
    }
}
// transform ray from local to world space
pub fn (f Frame) l2w_ray(r Ray) Ray {
    return Ray{
        e: f.l2w_point(r.e),
        d: f.l2w_direction(r.d),
        t_min: r.t_min,
        t_max: r.t_max,
    }
}

pub fn (v Vector) reflect(about Direction) Vector {
    return v.negate() + about.scale(2.0 * v.dot(about))
}
pub fn (d Direction) reflect(about Direction) Direction {
    return (d.negate().as_vector() + about.scale(2.0 * d.dot(about))).as_direction()
}



// Sphere movement helper
pub fn (f Frame) move_sphere(time f64) Frame {
    return Frame{
        o: f.o.add(f.velocity.scale(time)),
        x: direction_x
        y: direction_y
        z: direction_z

        velocity: f.velocity
    }
}




