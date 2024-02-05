# Circles and Lines

* [Demo](https://zotho.github.io/projects/circles_lines/)

[![Demo](https://raw.githubusercontent.com/zotho/zotho.github.io/master/projects/circles_lines/screenshot.png)](https://zotho.github.io/projects/circles_lines/)

* Control: mouse buttons

Written in Rust with [Macroquad](https://github.com/not-fl3/macroquad). Compiled to WASM

```console
# Build native
cargo run --example circles --release -- 100

# Build in WASM
cargo build \
    --example circles \
    --release \
    --target wasm32-unknown-unknown
cp target/wasm32-unknown-unknown/release/examples/circles.wasm circles.wasm
wasm-strip circles.wasm
python -m http.server
```