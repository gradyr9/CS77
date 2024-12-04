module gfx

import math

// IMPORTANT: image coordinate system has
// - x increasing to the right, and
// - y increasing down


// RGB Image
@[noinit]
pub struct Image {
pub:
    size Size2i @[required]
mut:
    data [][]Color      // row-major ordering
}

// RGBA Image
@[noinit]
pub struct Image4 {
    size Size2i @[required]
mut:
    data [][]Color4     // row-major ordering
}


///////////////////////////////////////////////////////////
// image generation functions

pub fn Image.new(size Size2i) Image {
    mut image := Image{ size: size }
    image.clear()
    return image
}

pub fn Image4.new(size Size2i) Image4 {
    mut image := Image4{ size: size }
    image.clear()
    return image
}

pub fn render_image0(size Size2i) Image4 {
    w, h := size.width, size.height
    wh, hh := f64(w) / 2.0, f64(h) / 2.0
    mrad := math.min(wh, hh)
    radii := [
        mrad * (205.0 / 256.0),
        mrad * (200.0 / 256.0),
        mrad * (150.0 / 256.0),
        mrad * (145.0 / 256.0),
        mrad * (100.0 / 256.0),
    ]
    mut image := Image4.new(size)
    mut r, mut g, mut b, mut a := f64(0), f64(0), f64(0), f64(0)
    for y in 0..h {
        for x in 0..w {
            rad := math.sqrt(math.pow(x - wh, 2.0) + math.pow(y - hh, 2.0))
                 if rad > radii[0] { r,g,b,a = 0.0, 1.0, 0.0, 0.0 }
            else if rad > radii[1] { r,g,b,a = 0.0, 0.0, 0.0, 1.0 }
            else if rad > radii[2] { r,g,b,a = 1.0, 1.0, 1.0, 1.0 }
            else if rad > radii[3] { r,g,b,a = 0.0, 0.0, 0.0, 1.0 }
            else if rad > radii[4] { r,g,b,a = 0.2, 0.2, 0.8, 0.5 }
            else                   { r,g,b,a = 1.0, 1.0, 1.0, (x+y)%2 }
            image.set_xy(x, y, Color4{ r, g, b, a })
        }
    }
    return image
}

pub fn render_image1(size Size2i, alpha bool) Image4 {
    w, h := size.width, size.height
    hh := f64(h) / 2.0
    mut image := Image4.new(size)
    mut r, mut g, mut b, mut a := f64(0), f64(0), f64(0), f64(0)
    for y in 0..h {
        for x in 0..w {
            v := math.abs(math.sin((y - hh) / hh * math.pi / 2.0))
            checker := int(math.fmod(
                math.floor(x * 20.0 / f64(w)) + math.floor(math.asin((y - hh) / hh) * 20.0 / (math.pi / 2.0)),
                2.0
            ))
            if checker == 0 { r,g,b = 1.0, 0.1, 0.1 }
            else            { r,g,b = 0.1, 1.0, 0.1 }
            a = if alpha { math.pow(v, 2.0) } else { 1.0 }
            image.set_xy(x, y, Color4{ r, g, b, a })
        }
    }
    return image
}


///////////////////////////////////////////////////////////
// convenience getter methods

pub fn (image Image) width() int {
    return image.size.width
}
pub fn (image Image) height() int {
    return image.size.height
}

pub fn (image Image4) width() int {
    return image.size.width
}
pub fn (image Image4) height() int {
    return image.size.height
}


///////////////////////////////////////////////////////////
// initializing and clearing methods

pub fn (mut i Image) clear() {
    i.data = [][] Color {
        len: i.size.height,
        init: [] Color {
            len: i.size.width
        }
    }
}
pub fn (mut i Image4) clear() {
    i.data = [][] Color4 {
        len: i.size.height,
        init: [] Color4 {
            len: i.size.width
        }
    }
}


///////////////////////////////////////////////////////////
// pixel setter and getter methods

pub fn (i Image) get(p Point2i) Color {
    return i.data[p.y][p.x]
}
pub fn (i Image) get_xy(x int, y int) Color {
    return i.data[y][x]
}
pub fn (i Image) get_color4(p Point2i) Color4 {
    return i.data[p.y][p.x].as_color4()
}
pub fn (i Image) get_xy_color4(x int, y int) Color4 {
    return i.data[y][x].as_color4()
}

pub fn (i Image4) get(p Point2i) Color4 {
    return i.data[p.y][p.x]
}
pub fn (i Image4) get_xy(x int, y int) Color4 {
    return i.data[y][x]
}
pub fn (i Image4) get_color(p Point2i) Color {
    return i.data[p.y][p.x].as_color()
}
pub fn (i Image4) get_xy_color(x int, y int) Color {
    return i.data[y][x].as_color()
}

pub fn (mut i Image) set(p Point2i, c Color) {
    i.set_xy(p.x, p.y, c)
}
pub fn (mut i Image) set_xy(x int, y int, c Color) {
    if x < 0 || x >= i.size.width { return }
    if y < 0 || y >= i.size.height { return }
    i.data[y][x] = c
}

pub fn (mut i Image4) set(p Point2i, c Color4) {
    i.set_xy(p.x, p.y, c)
}
pub fn (mut i Image4) set_xy(x int, y int, c Color4) {
    if x < 0 || x >= i.size.width { return }
    if y < 0 || y >= i.size.height { return }
    i.data[y][x] = c
}

pub fn (mut i Image4) set_color(p Point2i, c Color) {
    i.set_xy(p.x, p.y, c.as_color4())
}
pub fn (mut i Image4) set_xy_color(x int, y int, c Color) {
    if x < 0 || x >= i.size.width { return }
    if y < 0 || y >= i.size.height { return }
    i.data[y][x] = c.as_color4()
}
