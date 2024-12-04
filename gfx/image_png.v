module gfx

import stbi


///////////////////////////////////////////////////////////
// PNG loading and saving functions


pub fn load_png(filename string) Image4 {
    png := stbi.load(filename) or { panic(err) }
    size := png.width * png.height
    mut pixels := []Color4_u8{ len:size, cap:size }
    unsafe {
        pixels.data = png.data
    }
    mut image := Image4.new( Size2i{ png.width, png.height } )
    for y in 0 .. png.height {
        for x in 0 .. png.width {
            pixel := pixels[y * png.width + x]
            image.set_xy(x, y, pixel.as_color4())
        }
    }
    return image
}

pub fn (image Image) save_png(filename string) {
    size := image.width() * image.height()
    mut pixels := []Color4_u8{ len:size, cap:size }
    for y in 0 .. image.height() {
        for x in 0 .. image.width() {
            pixels[y * image.width() + x] = image.get_xy(x, y).as_color4_u8(1.0)
        }
    }
    stbi.stbi_write_png(filename, image.width(), image.height(), 4, pixels.data, image.width() * 4) or { panic(err) }
}

pub fn (image Image4) save_png(filename string) {
    size := image.width() * image.height()
    mut pixels := []Color4_u8{ len:size, cap:size }
    for y in 0 .. image.height() {
        for x in 0 .. image.width() {
            pixels[y * image.width() + x] = image.get_xy(x, y).as_color4_u8(1.0)
        }
    }
    stbi.stbi_write_png(filename, image.width(), image.height(), 4, pixels.data, image.width() * 4) or { panic(err) }
}

