module gfx

import os
import strings

///////////////////////////////////////////////////////////
// NewPBM loading and saving functions


@[params]
pub struct Load_NetPBM_Params {
    // filename to read from (full relative path)
    filename string @[required]
    // if not empty, the filename to read mask in as alpha values
    filename_mask string
}

pub fn load_netpbm_image(filename string) Image {
    image4 := load_netpbm_image4(filename:filename)
    size := image4.size
    mut image := Image.new(size)
    for y in 0 .. size.height {
        for x in 0 .. size.width {
            image.set_xy(x, y, image4.get_xy(x, y).as_color())
        }
    }
    return image
}

pub fn load_netpbm_image4(params Load_NetPBM_Params) Image4 {
    data := os.read_bytes(params.filename) or { panic(err) }
    mut netpbm := Loader{ data: data }
    assert netpbm.get() == `P`
    ver_num := netpbm.get_digit_ascii()
    assert 1 <= ver_num && ver_num <= 7

    mut width := 0
    mut height := 0
    mut depth := 0
    mut max_val := 0
    mut tuple_type := ''
    mut is_binary := false
    mut is_rgb := false
    mut has_alpha := false

    if ver_num == 7 {
        // PAM Loader, special case
        for {
            netpbm.eat_whitespace()
            attrib := netpbm.get_word()
            netpbm.get() // eat space (' ') or newline ('\n')
            match attrib {
                'ENDHDR'   { break }
                'WIDTH'    { width      = netpbm.get_int_ascii() }
                'HEIGHT'   { height     = netpbm.get_int_ascii() }
                'DEPTH'    { depth      = netpbm.get_int_ascii() }
                'MAXVAL'   { max_val    = netpbm.get_int_ascii() }
                'TUPLTYPE' { tuple_type = netpbm.get_string() }
                else { assert false }
            }
        }
        assert width > 0 && height > 0 && max_val > 0 && depth > 0
        assert tuple_type in ['RGB', 'RGB_ALPHA']
        has_alpha = tuple_type == "RGB_ALPHA"
        is_binary = true
    } else {
        // PBM / PGM / PPM Loader
        is_binary = ver_num >= 4
        is_rgb = ver_num == 3 || ver_num == 6
        depth = if is_rgb { 3 } else { 1 }
        has_alpha = false
        if ver_num % 3 == 1 { tuple_type = "BLACKANDWHITE" }
        if ver_num % 3 == 2 { tuple_type = "GRAYSCALE" }
        if ver_num % 3 == 0 { tuple_type = "RGB" }
        // process header
        vals_to_load := if tuple_type == "BLACKANDWHITE" { 2 } else { 3 }
        mut vals := []int{}
        for vals.len < vals_to_load {
            if netpbm.next_is_whitespace() { netpbm.eat_whitespace() }
            if netpbm.next_is_comment()    { netpbm.get_string() }
            if netpbm.next_is_digit()      { vals << netpbm.get_int_ascii() }
        }
        // eat single whitespace that separates header from data
        netpbm.get()
        width   = vals[0]
        height  = vals[1]
        max_val = if vals_to_load == 3 { vals[2] } else { 1 }
    }

    // not implemented, yet
    assert !(is_binary && max_val > 255)
    assert !(is_binary && ver_num == 1)
    // make sure all settings are as expected
    assert !(depth == 1 && has_alpha)
    assert !(depth == 2 && !has_alpha)
    assert !(depth == 3 && has_alpha)
    assert !(depth == 4 && !has_alpha)
    assert 1 <= depth && depth <= 4

    // use read_image4 instead
    //assert has_alpha

    size := Size2i{width, height}
    mut image := Image4.new(size)
    for y in 0 .. height {
        for x in 0 .. width {
            mut vs := []f64{}
            for _ in 0 .. depth {
                if is_binary {
                    vs << f64(netpbm.get_int_binary(1)) / f64(max_val)
                } else {
                    netpbm.eat_whitespace()
                    vs << f64(netpbm.get_int_ascii()) / f64(max_val)
                }
            }
            if depth == 1 && !has_alpha {
                image.set_xy(x, y, Color4{ vs[0], vs[0], vs[0], 1 })
            } else if depth == 2 && has_alpha {
                image.set_xy(x, y, Color4{ vs[0], vs[0], vs[0], vs[1] })
            } else if depth == 3 && !has_alpha {
                image.set_xy(x, y, Color4{ vs[0], vs[1], vs[2], 1 })
            } else if depth == 4 && has_alpha {
                image.set_xy(x, y, Color4{ vs[0], vs[1], vs[2], vs[3] })
            }
        }
    }

    if params.filename_mask != "" {
        image_mask := load_netpbm_image4(filename:params.filename_mask)
        assert image.size.width == image_mask.size.width
        assert image.size.height == image_mask.size.height
        // overwrite any alpha values with
        for y in 0 .. height {
            for x in 0 .. width {
                ci := image.get_xy(x, y)
                cm := image_mask.get_xy(x, y)
                image.set_xy(x, y, Color4{ ci.r, ci.g, ci.b, cm.r })
            }
        }
    }

    return image
}

@[params]
pub struct Save_NetPBM_Params {
    // image/image4 to write to file
    image  ?Image
    image4 ?Image4

    // filename to write to (full relative path)
    filename string @[required]
    // if not empty, the filename to write alpha values into (P2/P5)
    filename_mask string

    // max rgb value (clip)
    max f64 = 1.0

    // true:  write values in ascii  (P1/P2/P3)
    // false: write values in binary (P4/P5/P6)
    ascii bool

    // true: write alpha channel in image (always creates a P7)
    with_alpha bool

    // color to use as background when not writing alpha
    background Color = black
}

pub fn save_netpbm_image(f Save_NetPBM_Params) {
    mut w, mut h := 0, 0
    mut c := fn (x int, y int) []u8 { return [] }
    mut ca := fn (x int, y int) u8 { return 255 }

    if img := f.image {
        w, h = img.size.width, img.size.height
        data := img.data
        c = fn [f, data] (x int, y int) []u8 {
            return data[y][x].rgb_u8(f.max)
        }
    } else if img := f.image4 {
        w, h = img.size.width, img.size.height
        data := img.data
        if f.with_alpha {
            c = fn [f, data] (x int, y int) []u8 {
                return data[y][x].rgba_u8(f.max)
            }
        } else if f.filename_mask != '' {
            c = fn [f, data] (x int, y int) []u8 {
                return data[y][x].rgb_u8(f.max)
            }
        } else {
            c = fn [f, data] (x int, y int) []u8 {
                color4 := data[y][x]
                return (f.background.lerp(color4.as_color(), color4.a)).rgb_u8(f.max)
            }
        }
        ca = fn [data] (x int, y int) u8 {
            return data[y][x].a_u8()
        }
    }

    assert w != 0 && h != 0

    if f.with_alpha {
        mut data := []u8{}
        // write header
        data << 'P7\n'.bytes()
        data << 'WIDTH $w\n'.bytes()
        data << 'HEIGHT $h\n'.bytes()
        data << 'DEPTH 4\n'.bytes()
        data << 'MAXVAL 255\n'.bytes()
        data << 'TUPLTYPE RGB_ALPHA\n'.bytes()
        data << 'ENDHDR\n'.bytes()
        for y in 0..h {
            for x in 0..w {
                data << c(x, y)
            }
        }
        os.write_file_array(f.filename, data) or { panic(err) }
        return
    }

    if f.ascii {
        mut data := strings.new_builder(3*w*h*4)
        // write header
        data.write_string('P3\n')
        data.write_string('$w $h\n')
        data.write_string('255\n')
        for y in 0..h {
            for x in 0..w {
                rgb := c(x, y)
                r, g, b := rgb[0], rgb[1], rgb[2]
                data.write_string('$r $g $b\n')
            }
        }
        os.write_file(f.filename, data.str()) or { panic(err) }

        if f.filename_mask != '' {
            data = strings.new_builder(w*h*4)
            // write header
            data.write_string('P2\n$w $h\n255\n')
            for y in 0..h {
                for x in 0..w {
                    a := ca(x, y)
                    data.write_string('$a\n')
                }
            }
            os.write_file(f.filename_mask, data.str()) or { panic(err) }
        }
    } else {
        mut data := []u8{}
        // write header
        data << 'P6\n'.bytes()
        data << '$w $h\n'.bytes()
        data << '255\n'.bytes()
        for y in 0..h {
            for x in 0..w {
                data << c(x, y)
            }
        }
        os.write_file_array(f.filename, data) or { panic(err) }

        if f.filename_mask != '' {
            data = []
            // write header
            data << 'P5\n'.bytes()
            data << '$w $h\n'.bytes()
            data << '255\n'.bytes()
            for y in 0..h {
                for x in 0..w {
                    data << ca(x, y)
                }
            }
            os.write_file_array(f.filename_mask, data) or { panic(err) }
        }
    }
}

