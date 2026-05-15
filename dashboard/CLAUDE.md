# Yleaf Dashboard — Claude Instructions

## App icon

Source icon: `src-tauri/icons/icon.png` (1024×1024, committed to repo).

The icon is explicitly embedded into the binary via `include_bytes!` in
`src-tauri/src/lib.rs` and set at window startup with `set_icon()`. This
is required for the icon to appear correctly in MobaXterm's taskbar
(`_NET_WM_ICON` is not set automatically on Linux for X11-forwarded windows).

To update the icon: replace `src-tauri/icons/icon.png` with the new source,
then regenerate all sizes:
```bash
npm run tauri icon src-tauri/icons/icon.png
```
Then commit all changed files in `src-tauri/icons/` and rebuild.

**Never regenerate icons with ImageMagick directly** — use `npm run tauri icon`
which produces correctly formatted ICO/ICNS/PNG files that Tauri can load.

## Running in development

```bash
cd /media/disk1/arwin/yleaf-dashboard && npm run tauri dev
```

## Building for production

```bash
cd /media/disk1/arwin/yleaf-dashboard && npm run tauri build
```

AppImage output: `src-tauri/target/release/bundle/appimage/Yleaf 4.0_0.1.0_amd64.AppImage`
