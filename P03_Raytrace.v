module main
import os
import math
import gfx

////////////////////////////////////////////////////////////////////////////////////////
// Comment out lines in array below to prevent re-rendering every scene.
// If you create a new scene file, add it to the list below.
//
// NOTE: **BEFORE** you submit your solution, uncomment all lines, so
//       your code will render all the scenes!

const (
    scene_filenames = [
        'P03_00_motion_blur',
        'P03_01_depth_of_field',
        'P03_02_depth_of_field_2',
        'P03_03_blurry_reflections',
        'P03_04_creative_artifact',
    ]
)


////////////////////////////////////////////////////////////////////////////////////////
// module aliasing to make code a little easier to read
// ex: replacing `gfx.Scene` with just `Scene`

type Point     = gfx.Point
type Vector    = gfx.Vector
type Direction = gfx.Direction
type Normal    = gfx.Normal
type Ray       = gfx.Ray
type Color     = gfx.Color
type Image     = gfx.Image

type Intersection = gfx.Intersection
type Surface      = gfx.Surface
type Scene        = gfx.Scene


////////////////////////////////////////////////////////////////////////////////////////
// functions to implement


fn intersect_ray_surface(surface Surface, ray Ray) Intersection {
    /*
        if surface's shape is a sphere
            if ray does not intersect sphere, return no intersection
            compute ray's t value(s) at intersection(s)
            if ray's t is not a valid (between ray's min and max), return no intersection
            return intersection information
            NOTE: a ray can intersect a sphere in 0, 1, or 2 different points!

        if surface's shape is a quad
            if ray does not intersect plane, return no intersection
            compute ray's t value at intersection with plane
            if ray's t is not a valid (between min and max), return no intersection
            if intersection is outside the quad, return no intersection
            return intersection information
    */

    if surface.shape == gfx.Shape.sphere {
        time := ray.time
        sphere := surface.frame.move_sphere(time)


        ec := sphere.o.vector_to(ray.e) 
        a := 1.0 // always 1
        b := 2.0 * ray.d.dot(ec)
        c := ec.length_squared() - surface.radius * surface.radius
        determinant := b * b - 4 * a * c

        if determinant < 0 {
            return gfx.no_intersection
        }


        t1 := (-b - math.sqrt(determinant)) / (2.0 * a)
        t2 := (-b + math.sqrt(determinant)) / (2.0 * a)

        mut t := 0.0

        if !ray.valid_t(t1) && !ray.valid_t(t2) { 
            return gfx.no_intersection
        } else if ray.valid_t(t1) && ray.valid_t(t2) {
            t = math.min(t1, t2)
        } else if !ray.valid_t(t1) && ray.valid_t(t2) {
            t = t2
        } else if ray.valid_t(t1) && !ray.valid_t(t2) {
            t = t1
        }
        
        pt := ray.at(t)
        n := sphere.o.direction_to(pt)
        theta := math.acos(n.z)
        phi := math.atan2(n.y, n.x)
        x_dir := Vector{ math.sin(phi), math.cos(phi), 0 }.as_direction()
        y_dir := Vector{ math.cos(theta) * math.cos(phi), math.cos(theta) * math.sin(phi), math.sin(theta) }.as_direction()

        frame := gfx.Frame { 
            o: pt,
            x: x_dir,
            y: y_dir,
            z: n,
        }

        return Intersection { 
            distance : t, 
            surface: surface,
            frame: frame 
        }
    } else if surface.shape == gfx.Shape.quad {

        n := surface.frame.z

        if ray.d.dot(n) != 0 {
            ec := ray.e.vector_to(surface.frame.o)
            t := (ec.dot(n)) / (ray.d.dot(n))

            if ray.valid_t(t) {
                pt := ray.at(t)

                if surface.frame.o.vector_to(pt).linf_norm() <= surface.radius {
                    return Intersection {
                        distance: t, 
                        surface: surface, 
                        frame: gfx.Frame {
                            o: pt, 
                            x: surface.frame.x, 
                            y: surface.frame.y, 
                            z: n
                        }
                    }
                }
            }       
        }
    } else if surface.shape == gfx.Shape.circle {
        n := surface.frame.z

        if ray.d.dot(n) != 0 {
            ec := ray.e.vector_to(surface.frame.o)
            t := (ec.dot(n)) / (ray.d.dot(n))
            if ray.valid_t(t) {
                pt := ray.at(t)

                if surface.frame.o.vector_to(pt).l2_norm() <= surface.radius {
                    return Intersection {
                        distance: t, 
                        surface: surface, 
                        frame: gfx.Frame {
                            o: pt, 
                            x: surface.frame.x, 
                            y: surface.frame.y, 
                            z: n
                        }
                    }
                }
            }
        }       
    }


    return gfx.no_intersection
}

// Determines if given ray intersects any surface in the scene.
// If ray does not intersect anything, null is returned.
// Otherwise, details of first intersection are returned as an `Intersection` struct.
fn intersect_ray_scene(scene Scene, ray Ray) Intersection {
    mut closest := gfx.no_intersection  // type is Intersection

    /*
        for each surface in surfaces
            continue if ray did not hit surface ( ex: inter.miss() )
            continue if new intersection is not closer than previous closest intersection
            set closest intersection to new intersection
    */

    // Loop through each surface in the scene
    for surface in scene.surfaces {
        // Compute the intersection of the ray with the current surface
        intersection := intersect_ray_surface(surface, ray)

        // Check if there is an intersection and if it is closer than the current closest
        if !intersection.miss() && intersection.distance < closest.distance {
            // Update closest intersection
            closest = intersection 
        }
    }

    return closest  // return closest intersection
}

// Computes irradiance (as Color) from scene along ray
fn irradiance(scene Scene, ray Ray, filename string, depth int) Color {
    if depth > 5 {
        return gfx.black 
    }

    mut accum := gfx.black

    /*
        get scene intersection
        if not hit, return scene's background intensity
        accumulate color starting with ambient
        foreach light
            compute light response    (L)
            compute light direction   (l)
            compute light visibility  (V)
            compute material response (BRDF*cos) and accumulate
        if material has reflections (lightness of kr > 0)
            create reflection ray
            accumulate reflected light (recursive call) scaled by material reflection
        return accumulated color
    */

    // Get scene intersection
    intersection := intersect_ray_scene(scene, ray)

    // If not hit, return scene's background intensity
    if intersection.miss() {
        return scene.background_color
    }

    // Get the material of the intersected surface
    material := intersection.surface.material

    // Accumulate color starting with ambient
    accum += scene.ambient_color.mult(material.kd)
  
    // foreach light
    for light in scene.lights {
        // compute light response (L) and direction (l)
        s := light.frame.o
        p := intersection.frame.o
        p_to_s := p.vector_to(s)
        l := p_to_s.normalize()
        li := light.kl.scale(1.0 / (p_to_s.length_squared()))

        e := scene.camera.frame.o
        v := p.vector_to(e).normalize()

        if filename == 'P02_10_creative' {
            accum += gfx.Color{0.0, 0.0, 0.2} 
        }

        // compute light visibility (V)
        shadow_ray := Ray{
            e: p,
            d: l.direction(),
            //t_min = 0 is cool speckled effect used for creative artifact
        }
        shadow_intersection := intersect_ray_scene(scene, shadow_ray)
        

        // If the light is blocked by another surface, skip this light
        if !shadow_intersection.miss() && shadow_intersection.distance < p_to_s.length(){
            continue
        }


        n := intersection.frame.z
        h := l.add(v).normalize()

        pd := material.kd
        ps := material.ks.scale(math.pow(math.max(0.0, n.dot(h)), material.n))
        
        // Accumulate the color from this light
        accum += li.scale(math.abs(n.dot(l))).mult((pd + ps))
    }

    
    if intersection.surface.material.kr != gfx.black {
        fuzz := intersection.surface.material.fuzz
        mut reflection_color := gfx.black
        for _ in 0 .. scene.camera.sensor.blur_samples {
            origin_point := intersection.frame.o.add(intersection.frame.x.scale((f64(gfx.int_in_range(-500, 500)) / 1000.0) * fuzz)).add(intersection.frame.y.scale((f64(gfx.int_in_range(-500, 500)) / 1000.0) * fuzz)).add(intersection.frame.z.scale((f64(gfx.int_in_range(-500, 500)) / 1000.0) * fuzz))
            viewdir := origin_point.direction_to(scene.camera.frame.o)
            raydir := viewdir.scale(-1) + (intersection.frame.z.scale(2 * intersection.frame.z.dot(viewdir)))
            reflection_color.add_in(irradiance(scene, Ray { e: intersection.frame.o, d: raydir.as_direction()}, filename, depth + 1).mult(intersection.surface.material.kr))
        }
        reflection_color.scale_in(1.0 / f64(scene.camera.sensor.blur_samples))
        accum.add_in(reflection_color)
    }


    // Added in blue light to give everything a blue tinge for the creative Artifact
    if filename == 'P02_10_circle' {
        accum += gfx.Color{0.0, 0.0, 0.2} 
    }
    return accum
}


// Computes image of scene using basic Whitted raytracer.

fn raytrace(scene Scene, filename string) Image {
    mut image := gfx.Image.new(scene.camera.sensor.resolution)
    num_samples := scene.camera.sensor.samples
    dof_samples := scene.camera.sensor.dof_samples
    sensor_distance := scene.camera.sensor.distance
    camera_frame := scene.camera.frame
    resolution := scene.camera.sensor.resolution
    aperture_size := scene.camera.aperture_size

    // If no antialiasing
    if num_samples <= 1 {
        for y in 0 .. resolution.height {
            for x in 0 .. resolution.width {
                mut accum := gfx.black

                v := (f64(resolution.height - y) + 0.5) / f64(resolution.height)
                u := (f64(x) + 0.5) / f64(resolution.width)
                q_x := camera_frame.x.scale((u - 0.5) * resolution.width)
                q_y := camera_frame.y.scale((v - 0.5) * resolution.height)
                q_z := camera_frame.z.scale(-sensor_distance * resolution.width)
                q := scene.camera.frame.o.add(q_x + q_y + q_z)

                if dof_samples == 1 {
                    viewray := Ray { e: scene.camera.frame.o, d: scene.camera.frame.o.direction_to(q) }

                    irrad := irradiance(scene, viewray, filename, 0)
                    image.set_xy(y, x, irrad)
                }
                else {
                    // With DOF
                    for _ in 0 .. dof_samples {
                        base_direction := scene.camera.frame.o.direction_to(q)
                        focal_point := camera_frame.o.add(base_direction.scale(scene.camera.focal_distance))
                        lens_point := camera_frame.o.add(camera_frame.x.scale((f64(gfx.int_in_range(-5, 5)) / 10.0) * aperture_size)).add(camera_frame.y.scale((f64(gfx.int_in_range(-5, 5)) / 10.0) * aperture_size))

                        viewray := Ray { 
                            e: lens_point, 
                            d: lens_point.direction_to(focal_point) 
                        }

                        irrad := irradiance(scene, viewray, filename, 0)
                        accum = accum.add(irrad)
                    }
                    accum.scale_in(1.0 / f64(dof_samples))
                    image.set_xy(x, y, accum)
                }
            }
        }
    } else {
        // If antialiasing and depth of field
        for y in 0 .. resolution.height {
            for x in 0 .. resolution.width {
                mut accum := gfx.black

                for col in 0 .. num_samples {
                    for row in 0 .. num_samples {
                        // Jittered sampling for anti-aliasing
                        v := (f64(resolution.height - y) + f64(col)/f64(num_samples)) / f64(resolution.height)
                        u := (f64(x) + f64(row)/f64(num_samples)) / f64(resolution.width)
                        
                        qx := camera_frame.x.scale((u - 0.5) * resolution.width)
                        qy := camera_frame.y.scale((v - 0.5) * resolution.height)
                        qz := camera_frame.z.scale((-1) * sensor_distance * resolution.width)

                        q := scene.camera.frame.o.add(qx + qy + qz)

                        rand_time := f64(gfx.int_in_range(0, 1000)) / 1000.0

                        // Perform lens sampling for DoF
                        mut dof_accum := gfx.black
                        for _ in 0 .. dof_samples {
                            // Calculate focal point based on base direction
                            base_direction := camera_frame.o.direction_to(q)
                            focal_point := camera_frame.o.add(base_direction.scale(scene.camera.focal_distance))
                            
                            // Sample a random point on the lens aperture
                            lens_point := camera_frame.o.add(camera_frame.x.scale((f64(gfx.int_in_range(-500, 500)) / 1000.0) * aperture_size)).add(camera_frame.y.scale((f64(gfx.int_in_range(-500, 500)) / 1000.0) * aperture_size))
                            
                            
                            
                            // Create a new ray with DoF effect
                            viewray := Ray { 
                                e: lens_point, 
                                d: lens_point.direction_to(focal_point),
                                time: rand_time
                            }
                            
                            // Calculate irradiance for this DoF sample
                            irrad := irradiance(scene, viewray, filename, 0)
                            dof_accum = dof_accum.add(irrad)
                        }
                        
                        // Average the results of the DoF samples
                        dof_accum.scale_in(1.0 / f64(dof_samples))
                        
                        // Accumulate the result for anti-aliasing
                        accum = accum.add(dof_accum)
                    }
                }

                // Average the results of the anti-aliasing samples
                accum.scale_in(1.0 / math.pow(scene.camera.sensor.samples, 2))
                image.set_xy(x, y, accum)
            }
        }
    }
    return image
}





fn main() {
    // Make sure images folder exists, because this is where all generated images will be saved
    if !os.exists('output') {
        os.mkdir('output') or { panic(err) }
    }

    for filename in scene_filenames {
        println('Rendering ${filename}...')
        scene := gfx.scene_from_file('scenes/${filename}.json')!
        image := raytrace(scene, filename)
        image.save_png('output/${filename}.png')
    }

    println('Done!')
}
