module gfx

import math

pub struct Intersection {
pub:
    frame    Frame              // coordinate frame of intersection
                                // note: frame.o is the point of intersection, and
                                //       frame.z is aligned with normal at point of intersection

    surface  Surface            // surface that was hit.  should be an optional type, but that destroys
                                // readability of implementation, so use finite distance to indicate a hit

    distance f64 = math.inf(1)  // distance from ray origin to point of intersection
                                // note: infinite distance means no hit; finite distance means hit!
}


///////////////////////////////////////////////////////////
// a no-hit intersection const

pub const (
    no_intersection = Intersection{}
)


///////////////////////////////////////////////////////////
// convenience getters

pub fn (inter Intersection) o() Point {
    return inter.frame.o
}
pub fn (inter Intersection) normal() Normal {
    return inter.frame.z.as_normal()
}
pub fn (inter Intersection) material() Material {
    return inter.surface.material
}


///////////////////////////////////////////////////////////
// testing and comparison methods

// true if intersection distance is finite (hit)
pub fn (inter Intersection) hit() bool {
    return math.is_finite(inter.distance)
}

// true if intersection distance is infinitely far away (miss)
pub fn (inter Intersection) miss() bool {
    return math.is_inf(inter.distance, 1)
}

// returns true if `a` is closer than `b`
pub fn (a Intersection) is_closer(b Intersection) bool {
    return a.distance < b.distance
}
// found bug in V on 2022.08.15
// when using module import aliasing, the comparator operators cause problems
// pub fn (a Intersection) < (b Intersection) bool {
//     return a.distance < b.distance
// }
// pub fn (a Intersection) == (b Intersection) bool {
//     return a.distance == b.distance
// }


