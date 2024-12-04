module gfx

pub struct Vertex {
pub:
    point Point2i @[required]
    color Color
}

pub struct GeoPrimitive {
pub:
    p0 Point2i
    c0 Color
    p1 Point2i
    c1 Color
}

pub struct Fragment {
pub:
    point Point2i @[required]
    color Color   @[required]
}






