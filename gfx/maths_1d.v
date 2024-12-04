module gfx

import math
import rand

pub fn min3<T>(a T, b T, c T) T {
    return math.min(a, math.min(b, c))
}

pub fn max3<T>(a T, b T, c T) T {
    return math.max(a, math.max(b, c))
}

pub fn int_in_range(min int, max int) int {
    return rand.int_in_range(
        math.min(min, max),
        math.max(min, max),
    ) or { return min }
}

pub fn sign(v int) int {
    if v < 0 { return -1 }
    if v > 0 { return 1 }
    return 0
}


pub struct InclusiveRange {
    start int
    end int
mut:
    initialized bool
    current int
    delta int
}

pub fn (mut ir InclusiveRange) next() !int {
    if !ir.initialized {
        ir.initialized = true
        ir.current = ir.start
        ir.delta = sign(ir.end - ir.start)
    } else if ir.current == ir.end {
        return error('done iterating')
    } else {
        ir.current += ir.delta
    }
    return ir.current
}

