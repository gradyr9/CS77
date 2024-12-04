module gfx

import strings

struct Loader {
    data []u8
mut:
    offset int
}

fn (mut loader Loader) peek() u8 {
    return loader.data[loader.offset]
}
fn (mut loader Loader) get() u8 {
    d := loader.peek()
    loader.offset += 1
    return d
}

fn (mut loader Loader) next_is_digit() bool {
    match loader.peek() {
        `0` ... `9` { return true }
        else { return false }
    }
}
fn (mut loader Loader) next_is_char() bool {
    match loader.peek() {
        `a` ... `z` { return true }
        `A` ... `Z` { return true }
        else { return false }
    }
}
fn (mut loader Loader) next_is_comment() bool {
    return loader.peek() == `#`
}
fn (mut loader Loader) next_is_whitespace() bool {
    match loader.peek() {
        ` `  { return true }
        `\n` { return true }
        `\t` { return true }
        `\r` { return true }
        else { return false }
    }
}
fn (mut loader Loader) next_is_newline() bool {
    return loader.peek() == `\n`
}

fn (mut loader Loader) get_digit_ascii() int {
    return int(loader.get() - `0`)
}

fn (mut loader Loader) get_int_ascii() int {
    mut v := 0
    for loader.next_is_digit() {
        v = v * 10 + loader.get_digit_ascii()
    }
    return v
}
fn (mut loader Loader) get_int_binary(bytes int) int {
    mut v := 0
    for _ in 0..bytes {
        v = v * 256 + int(loader.get())
    }
    return v
}

fn (mut loader Loader) eat_whitespace() {
    for loader.next_is_whitespace() {
        loader.get()
    }
}

fn (mut loader Loader) get_word() string {
    mut s := strings.new_builder(10)
    for !loader.next_is_whitespace() {
        s.write_u8(loader.get())
    }
    return s.str()
}
fn (mut loader Loader) get_string() string {
    mut s := strings.new_builder(100)
    for !loader.next_is_newline() {
        s.write_u8(loader.get())
    }
    return s.str()
}
