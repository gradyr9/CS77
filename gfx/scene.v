module gfx

import os
import json


///////////////////////////////////////////////////////////
// structs used to define a scene
//
// NOTE: if any post-load processing needs to happen on
//       structs after loading from JSON file, add
//       processing code to `update_after_load` below!
//
//       ex: Frames can have Look At parameters (eye, target, up)
//           specified instead of usual Frame params (o, x, y, z).

pub struct Scene {
pub:
    camera           Camera
    background_color Color       = Color{ 0.2, 0.2, 0.2 }
    ambient_color    Color       = Color{ 0.2, 0.2, 0.2 }
    lights           [] Light    = [ Light{} ]
    surfaces         [] Surface  = [ Surface{} ]
}

pub struct Camera {
pub:
    sensor Sensor = Sensor{}
    frame  Frame  = Frame{
        o: Point{ 0.0, 0.0, 1.0 }
        x: direction_x
        y: direction_y
        z: direction_z
    }
    aperture_size f64 = 0.05  // Size of the aperture (controls blur)
    focal_distance f64 = 5.0 // Distance at which objects are in focus
}


pub struct Sensor {
pub:
    size       Size2  = Size2{  1.0, 1.0 }
    resolution Size2i = Size2i{ 512, 512 }
    distance   f64    = 1.0
    samples    int    = 1
    dof_samples int    = 1
    blur_samples int    = 1
}

pub struct Light {
pub:
    kl     Color = white
    frame  Frame = Frame{
        o: Point{ 0.0, 0.0, 5.0 }
        x: direction_x
        y: direction_y
        z: direction_z
    }
}

pub enum Shape {
    sphere
    quad
    circle
}

pub struct Surface {
pub:
    shape    Shape = Shape.sphere
    radius   f64   = 1.0
    frame    Frame
    material Material
}

pub struct Material {
pub:
    kd Color = white
    ks Color = black
    n  f64   = 10.0
    kr Color = black
    fuzz f64 = 0.0
}


// function to update scene after loading from JSON
fn update_after_load(scene Scene) Scene {
    // NOTE: attempted a more generic version using reflection (see bottom)
    return Scene{
        ...scene                                            // copy everything from scene with some overrides (below)
        camera: Camera{
            ...scene.camera                                 // copy everything from camera with override
            frame: scene.camera.frame.as_frame()            // update frame for camera if look at is specified
        }
        lights: scene.lights.map(fn (light Light) Light {   // for each light...
            return Light{
                ...light                                    // copy everything from camera with override
                frame: light.frame.as_frame()               // update frame for light if look at is specified
            }
        })

    }
}



///////////////////////////////////////////////////////////
// scene importing and exporting functions

pub fn scene_from_file(path string) !Scene {
    // load and decode scene from JSON file
    data := os.read_file(path)!
    scene := json.decode(Scene, data)!
    return update_after_load(scene)
}

pub fn (scene Scene) to_json() string {
    return json.encode(scene)
}



///////////////////////////////////////////////////////////
// convenience getters

pub fn (light Light) o() Point {
    return light.frame.o
}

pub fn (surface Surface) o() Point {
    return surface.frame.o
}

