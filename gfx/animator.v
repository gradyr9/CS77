module gfx

import gg
import gx
import os
import json
import time

const help_doc = '
usage:
    animator <animator_data.json>
  or
    animator [--title "<window title>"] [--size <width> <height>]
             [--background <r> <g> <b>] [--framerate <fps>]
             [--cycle <wrap|bounce>] [<png files>]
'

enum Cycle {
    wrap
    bounce
}

struct AnimatorData {
mut:
    title      string                                   // window title
    size       []int    = [ 512, 512 ]                  // window size
    background []u8     = [ u8(255), u8(0), u8(255) ]   // window background color (for transparent colors)
    framerate  int      = 30                            // framerate (fps) of animation
    cycle      Cycle    = .wrap                         // cycling option
    filenames  []string = []                            // filenames of PNG images

    // internal, not part of json
    start      time.Time    // timestamp when program started
    image_inds []int        // cache indices for images
}

struct App {
mut:
    start time.Time
    gg    &gg.Context = unsafe { nil }

    anim AnimatorData
    framerate  int
    wrap       bool
    image_inds []int
}

fn get_animator_data() !AnimatorData {
    if os.args.len < 2 {
        println(help_doc)
        exit(0)
    }

    mut anim := AnimatorData{}

    mut i := 1
    for i < os.args.len {
        arg := os.args[i]
        match arg {
            "--title" {
                anim.title = os.args[i+1]
                i += 2
            }
            "--size" {
                anim.size = [ os.args[i+1].int(),  os.args[i+2].int() ]
                i += 3
            }
            "--background" {
                anim.background = [ os.args[i+1].u8(), os.args[i+2].u8(), os.args[i+3].u8() ]
                i += 4
            }
            "--framerate" {
                anim.framerate = os.args[i+1].int()
                i += 2
            }
            "--cycle" {
                match os.args[i+1] {
                    "wrap" {
                        anim.cycle = .wrap
                    }
                    "bounce" {
                        anim.cycle = .bounce
                    }
                    else {
                        eprintln('Error: --cycle must be wrap or bounce')
                        exit(1)
                    }
                }
                i += 2
            }
            else {
                match os.file_ext(arg) {
                    '.json' {
                        data := os.read_file(arg)!
                        anim = json.decode(AnimatorData, data)!
                        i += 1
                    }
                    '.png' {
                        anim.filenames << arg
                        i += 1
                    }
                    else {
                        eprintln('Error: unexpected arg ${arg}')
                        exit(1)
                    }
                }
            }
        }
    }

    // validate correct arguments
    if anim.filenames.len == 0 {
        eprintln('Error: no PNGs specified')
        exit(1)
    }
    if anim.framerate <= 0 {
        eprintln('Error: invalid framerate (${anim.framerate} <= 0)')
        exit(1)
    }
    if anim.background.len != 3 {
        eprintln('Error: invalid background color (${anim.background}.len != 3)')
        exit(1)
    }
    if anim.size.len != 2 || anim.size[0] <= 0 || anim.size[1] <= 0 {
        eprintln('Error: invalid window size (${anim.size})')
        exit(1)
    }

    return anim
}

fn main() {
    mut anim := get_animator_data() or {panic('no animator data')}

    mut app := &App{}
    app.gg = gg.new_context(
        bg_color: gx.rgb(anim.background[0], anim.background[1], anim.background[2])
        width: anim.size[0]
        height: anim.size[1]
        create_window: true
        window_title: if anim.title == '' { 'Animator' } else { 'Animator: ${anim.title}' }
        frame_fn: frame
        user_data: app
    )

    // must create gg Context _before_ we start loading images!
    for path in anim.filenames {
        image := app.gg.create_image(path) or {
            println('Error: cannot load image ${path}')
            exit(1)
        }
        anim.image_inds << app.gg.cache_image(image)
    }

    anim.start = time.now()

    app.anim = anim
    app.gg.run()
}

fn get_index(anim AnimatorData) int {
    delta := (time.now() - anim.start).seconds()
    frame := int(anim.framerate * delta)
    indices := anim.image_inds.len

    match anim.cycle {
        .wrap {
            return frame % indices
        }
        .bounce {
            mut index := frame % (2 * indices - 2)
            if index >= indices {
                index = (indices - 2) - (index - indices)
            }
            return index
        }
    }
}

fn get_gg_image_index(anim AnimatorData) int {
    return anim.image_inds[get_index(anim)]
}

fn frame(mut app App) {
    index := get_gg_image_index(app.anim)
    image := app.gg.get_cached_image_by_idx(index)

    app.gg.begin()
    app.gg.draw_image(0, 0, image.width, image.height, image)
    app.gg.end()
}
