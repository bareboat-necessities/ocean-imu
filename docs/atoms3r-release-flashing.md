# Flash release binaries to an AtomS3R from the command line

This guide installs a prebuilt `ocean-imu` sketch on an **M5Stack AtomS3R
(ESP32-S3, 8 MB flash)** over USB. No Arduino IDE, repository checkout, sketch
compilation, or M5Unified/Eigen library installation is needed.

Use **esptool with the merged firmware** for the simplest installation. The
[Arduino CLI alternative](#alternative-upload-with-arduino-cli) uploads the
separate images and also covers archives without a merged image. Choose one
upload method, not both.

> **Back up before flashing.** Uploading replaces the running firmware. The
> merged raw image also writes the padding between its components, which can
> erase NVS and saved IMU calibration even without an `erase-flash` command.
> Plan to recalibrate afterward. To retain calibration during an update with
> the **same partition layout**, use the Arduino CLI alternative instead and
> still make a backup. Do not use `erase-flash` for a routine update.

## 1. Choose the release and sketch

Open the project's [GitHub Releases](https://github.com/bareboat-necessities/ocean-imu/releases).
The rolling [vTest prerelease](https://github.com/bareboat-necessities/ocean-imu/releases/tag/vTest)
receives successful `main` build uploads; a numbered release is a separate
choice. A published binary is not necessarily from the current `main` commit.

Expand **Assets** and choose the ZIP for exactly one sketch:

| Sketch name (`SKETCH` below) | Application |
| --- | --- |
| `atomS3R_imu_m5_basic` | Raw/calibrated IMU diagnostics |
| `atomS3R_mag_diagnostics` | Guided magnetometer/compass diagnostics |
| `atomS3R_compass_mahony` | Mahony compass / AHRS |
| `atomS3R_compass_qmekf` | Quaternion-MEKF compass / AHRS |
| `atomS3R_ins_kalman_ou2` | OU-II marine INS |
| `atomS3R_ins_kalman_ou3` | OU-III marine INS |
| `atomS3R_ins_tfg` | Two-frame Lie-group marine INS |
| `atomS3R_ins_pii_observer` | PII marine motion observer |
| `atomS3R_ins_nlo` | Time-varying-gain nonlinear observer |

Asset names have the form `atomS3R_ins_kalman_ou3_bin-YYYY-MM-DD.zip`.
Use the **actual filename in Assets**, not today's date. The rolling release
can retain several dates: select the newest matching ZIP unless deliberately
installing an older build. Do not select GitHub's **Source code (zip)** or a PDF.

The commands below use OU-III as an example. Change `SKETCH`/`$Sketch` for another
row and `RELEASE`/`$Release` for another release tag. Keep using the same terminal
and working directory through the following steps; stop if any command fails.

## 2. Install Python and esptool

[esptool 5.x](https://docs.espressif.com/projects/esptool/en/latest/esp32s3/installation.html)
requires Python 3.10 or newer. The version range below keeps the commands on
esptool 5.x, whose command names use hyphens (`write-flash`, `read-flash`).
A virtual environment avoids changing system Python or Arduino's bundled tools.

### Windows PowerShell

With [WinGet](https://learn.microsoft.com/en-us/windows/dev-environment/python)
available, install Python if needed:

```powershell
winget install --exact --id Python.Python.3.14
```

Close and reopen PowerShell after installation. If WinGet is unavailable, use
[Python's Windows installer](https://www.python.org/downloads/windows/).
Then run:

```powershell
python --version
New-Item -ItemType Directory -Force "$HOME\ocean-imu-flash" | Out-Null
Set-Location "$HOME\ocean-imu-flash"
python -m venv .venv
$Python = Join-Path $PWD '.venv\Scripts\python.exe'
& $Python -m pip install --upgrade pip
& $Python -m pip install 'esptool>=5,<6'
& $Python -m esptool version
```

Use `py -3` instead of `python` for the first two Python commands if that is how
Python is installed on your machine. No virtual-environment activation or
PowerShell execution-policy change is necessary.

### macOS / Linux (Bash)

On macOS, install Python with [Homebrew](https://brew.sh/) if needed:

```bash
brew install python
```

On Debian/Ubuntu Linux, install the prerequisites if needed:

```bash
sudo apt-get update
sudo apt-get install -y python3 python3-venv python3-pip curl unzip
```

Then, on either system, use Bash and run:

```bash
python3 --version
mkdir -p "$HOME/ocean-imu-flash"
cd "$HOME/ocean-imu-flash"
python3 -m venv .venv
PYTHON="$PWD/.venv/bin/python"
"$PYTHON" -m pip install --upgrade pip
"$PYTHON" -m pip install 'esptool>=5,<6'
"$PYTHON" -m esptool version
```

Other Linux distributions need equivalent packages. Check that `python3` is
3.10 or newer; do not use `sudo pip` or bypass an externally-managed-Python warning.

## 3. Download and unzip the selected asset

Paste the ZIP's **filename**, not its URL, at the prompt. Use a new extraction
directory for each download so files from different sketches/builds cannot mix.
The commands intentionally do not overwrite an existing extraction directory.

### Windows PowerShell

```powershell
$ErrorActionPreference = 'Stop'
$Release = 'vTest'
$Sketch = 'atomS3R_ins_kalman_ou3'
$Asset = Read-Host 'Paste the exact ZIP asset filename from GitHub Releases'
if ($Asset -notlike "${Sketch}_bin-*.zip") { throw "Select the ZIP for $Sketch." }
$Url = "https://github.com/bareboat-necessities/ocean-imu/releases/download/$Release/$Asset"
$BinDir = Join-Path $PWD ([IO.Path]::GetFileNameWithoutExtension($Asset))
if (Test-Path $BinDir) { throw 'Use a new extraction directory for this download.' }
Invoke-WebRequest -Uri $Url -OutFile $Asset
Expand-Archive -LiteralPath $Asset -DestinationPath $BinDir
Get-ChildItem -LiteralPath $BinDir
```

### macOS / Linux (Bash)

```bash
RELEASE='vTest'
SKETCH='atomS3R_ins_kalman_ou3'
read -r -p 'Paste the exact ZIP asset filename from GitHub Releases: ' ASSET
[[ "$ASSET" == "${SKETCH}_bin-"*.zip ]] || { echo "Select the ZIP for $SKETCH."; exit 1; }
URL="https://github.com/bareboat-necessities/ocean-imu/releases/download/$RELEASE/$ASSET"
BIN_DIR="$PWD/${ASSET%.zip}"
[[ ! -e "$BIN_DIR" ]] || { echo 'Use a new extraction directory for this download.'; exit 1; }
curl --fail --location --retry 3 "$URL" -o "$ASSET" &&
  unzip -t "$ASSET" &&
  mkdir "$BIN_DIR" &&
  unzip "$ASSET" -d "$BIN_DIR"
ls -l "$BIN_DIR"
```

The current [release packaging workflow](../.github/workflows/build.yml) puts
these files directly in the ZIP root (names below use the chosen sketch):

| File | Use |
| --- | --- |
| `<sketch>_firmware.bin` | Merged image: flash at **`0x0`** with esptool |
| `<sketch>.ino.bin` | Application only: **not** a complete image for `0x0` |
| `<sketch>.ino.bootloader.bin` | Separate ESP32-S3 bootloader |
| `<sketch>.ino.partitions.bin` | Compiled partition table |
| Other `.bin` / `.csv` files, when present | Build outputs; not extra files to flash blindly |

The workflow's merged image contains the bootloader at `0x0`, partition table at
`0x8000`, `boot_app0.bin` at `0xe000`, and application at `0x10000`, with DIO,
80 MHz flash, and an 8 MB flash-size header. It already includes `boot_app0.bin`;
that file does not have to be present separately in the ZIP. Do **not** rename
an application-only `.ino.bin` to `_firmware.bin`, flash every `.bin` with a
wildcard, or mix components from different builds.

## 4. Connect the AtomS3R and find its port

Use a USB-C **data** cable. Close Arduino Serial Monitor, terminal programs,
M5Burner, and any other program using the device's serial port.

Follow M5Stack's [AtomS3R download-mode instructions](https://docs.m5stack.com/en/core/AtomS3R#download-mode):
with USB connected, hold the **reset button** for about two seconds until the
internal green LED lights, then release it. The green LED goes out in download
mode. Use the reset button, not the display's application button.

List ports **after** entering download mode:

```powershell
# Windows PowerShell
& $Python -m serial.tools.list_ports -v
$Port = 'COM5'   # Replace with this AtomS3R's actual port.
& $Python -m esptool --chip esp32s3 --port $Port flash-id
```

```bash
# macOS / Linux (Bash)
"$PYTHON" -m serial.tools.list_ports -v
PORT='/dev/ttyACM0'   # Linux example; replace with the actual port.
# On macOS this is typically /dev/cu.usbmodem... instead.
"$PYTHON" -m esptool --chip esp32s3 --port "$PORT" flash-id
```

Disconnect/reconnect the board and compare the lists if several ports are
present. Confirm that esptool identifies an ESP32-S3 and 8 MB flash. Do not
force a chip mismatch. USB port names can change between download mode and the
running application; recheck the port whenever reconnecting or resetting.

### Back up the existing flash

Before overwriting a configured device, save its entire flash, including NVS.
Choose a new backup filename each time and keep it outside the extraction folder.
These commands read the device; they do not flash it:

```powershell
$Backup = 'atoms3r-backup-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '.bin'
& $Python -m esptool --chip esp32s3 --port $Port --baud 460800 read-flash 0 ALL $Backup
```

```bash
BACKUP="atoms3r-backup-$(date +%Y%m%d-%H%M%S).bin"
"$PYTHON" -m esptool --chip esp32s3 --port "$PORT" --baud 460800 read-flash 0 ALL "$BACKUP"
```

Check that the read succeeded; an 8 MB full-flash backup is 8,388,608 bytes.
A whole-device backup is for recovery of that same device, not a calibration
file to copy into another board or a different partition layout. Re-enter
download mode and check the port again if the tool reset the device.

## 5. Upload the merged firmware with esptool

**This method can erase saved calibration; see the warning at the top.**
Check for the exact `_firmware.bin` filename, then upload that one file at `0x0`.
The command keeps the flash settings already encoded by the release workflow.

### Windows PowerShell

```powershell
$Firmware = Join-Path $BinDir "${Sketch}_firmware.bin"
if (-not (Test-Path -LiteralPath $Firmware)) { throw 'No merged image: use the Arduino CLI alternative.' }
& $Python -m esptool --chip esp32s3 --port $Port --baud 460800 write-flash 0x0 $Firmware
```

### macOS / Linux (Bash)

```bash
FIRMWARE="$BIN_DIR/${SKETCH}_firmware.bin"
[[ -f "$FIRMWARE" ]] || { echo 'No merged image: use the Arduino CLI alternative.'; exit 1; }
"$PYTHON" -m esptool --chip esp32s3 --port "$PORT" --baud 460800 write-flash 0x0 "$FIRMWARE"
```

Let the upload and verification finish without unplugging the board. esptool
verifies the data it writes. After a successful upload, briefly press reset or
reconnect USB if the application does not start automatically.

## Alternative: upload with Arduino CLI

This uploads the **separate release images**, not the merged `_firmware.bin`.
It does not compile source. The Arduino ESP32 core supplies its own upload tool
and `boot_app0.bin`, so Arduino IDE, ESP-IDF, and sketch libraries are unnecessary.

For archives produced by the current workflow, install **`esp32:esp32@3.3.7`**
and use exactly this board configuration:

```text
esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc
```

The `m5stack_atoms3` identifier is intentional: it is the target used to produce
these AtomS3R release binaries. For a historical release, check its version of
[`.github/workflows/build.yml`](../.github/workflows/build.yml) and match the
core/board settings from that release rather than assuming today's settings.
Installing a different board profile does not change a precompiled binary.

### Install Arduino CLI

On **Windows x64 PowerShell**, from the same working directory:

```powershell
Invoke-WebRequest -Uri 'https://downloads.arduino.cc/arduino-cli/arduino-cli_latest_Windows_64bit.zip' -OutFile 'arduino-cli.zip'
Expand-Archive -LiteralPath 'arduino-cli.zip' -DestinationPath '.\arduino-cli-tools' -Force
$env:Path = (Join-Path $PWD 'arduino-cli-tools') + ';' + $env:Path
arduino-cli version
```

On **macOS / Linux**, install through Arduino's script (or use
`brew install arduino-cli` when Homebrew is already installed):

```bash
curl --fail --location https://raw.githubusercontent.com/arduino/arduino-cli/master/install.sh -o install-arduino-cli.sh
# Review the official installer before running it.
sh install-arduino-cli.sh
export PATH="$PWD/bin:$PATH"
arduino-cli version
```

The PATH changes above last for this terminal session. For other CPU/OS variants,
use the matching package from [Arduino CLI installation](https://docs.arduino.cc/arduino-cli/installation/).

### Install the matching ESP32 core and upload

These commands work in **both PowerShell and Bash**; the explicit additional
URL avoids changing an existing Arduino configuration:

```text
arduino-cli core update-index --additional-urls https://raw.githubusercontent.com/espressif/arduino-esp32/gh-pages/package_esp32_index.json
arduino-cli core install esp32:esp32@3.3.7 --additional-urls https://raw.githubusercontent.com/espressif/arduino-esp32/gh-pages/package_esp32_index.json
arduino-cli core list
arduino-cli board list
```

Confirm the board is in download mode and update `$Port`/`PORT` from the port
list. The board may be listed as unknown; the explicit FQBN supplies its profile.
Keep `<sketch>.ino.bin`, `<sketch>.ino.bootloader.bin`, and
`<sketch>.ino.partitions.bin` together in the extraction directory, with their
original names. Stop if any is missing; download a complete matching archive.

**Windows PowerShell:**

```powershell
$App = Join-Path $BinDir "$Sketch.ino.bin"
foreach ($Suffix in @('.ino.bin', '.ino.bootloader.bin', '.ino.partitions.bin')) {
    if (-not (Test-Path -LiteralPath (Join-Path $BinDir "$Sketch$Suffix"))) { throw "Missing $Sketch$Suffix" }
}
arduino-cli upload --fqbn 'esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc' --port $Port --input-file $App --upload-property 'upload.speed=460800' --verbose
```

**macOS / Linux (Bash):**

```bash
for suffix in .ino.bin .ino.bootloader.bin .ino.partitions.bin; do
  [[ -f "$BIN_DIR/$SKETCH$suffix" ]] || { echo "Missing $SKETCH$suffix"; exit 1; }
done
arduino-cli upload --fqbn 'esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc' --port "$PORT" --input-file "$BIN_DIR/$SKETCH.ino.bin" --upload-property 'upload.speed=460800' --verbose
```

`--input-file` identifies the application and its companion files unambiguously;
do not point it at `_firmware.bin` (Arduino's application offset is `0x10000`,
not `0x0`). The [ESP32 upload recipe](https://github.com/espressif/arduino-esp32/blob/3.3.7/platform.txt)
writes bootloader, partition table, `boot_app0.bin`, and application separately.
With an unchanged partition layout and no full erase, this avoids overwriting
the NVS gap as the merged raw image does. It is not a guarantee of calibration
format compatibility across releases; retain the backup and recalibrate when
required by the application.

## 6. Check the running sketch

After reset, find the application's USB serial port again. The compass/INS
sketches have an on-device UI; the basic IMU sketch prints diagnostics. Follow
the chosen sketch's README and any calibration prompts before evaluating it.
To view serial output at 115200 baud, use either installed tool:

```powershell
& $Python -m serial.tools.miniterm $Port 115200
# Or, after installing Arduino CLI:
arduino-cli monitor --port $Port --config baudrate=115200
```

```bash
"$PYTHON" -m serial.tools.miniterm "$PORT" 115200
# Or, after installing Arduino CLI:
arduino-cli monitor --port "$PORT" --config baudrate=115200
```

Use **Ctrl+]** to exit miniterm, or **Ctrl+C** to exit Arduino CLI monitor.
Close the monitor before another upload. Prebuilt firmware uses its compiled
settings; changing NMEA/UI options or other compile-time switches requires a
source build, as described in the [sensor Arduino setup](../sensors/README.md#arduino-installation).

## Troubleshooting

- **404 / wrong ZIP:** copy the exact tag and dated asset filename from Releases.
  Do not invent the date or use GitHub's source archive. An older archive without
  `_firmware.bin` needs the separate-image Arduino CLI method.
- **No port / failed connection:** try a known data cable and direct USB port,
  close serial programs, enter download mode with the reset button, and list
  ports again. Use native Windows PowerShell rather than WSL unless USB has
  explicitly been forwarded into WSL.
- **Linux permission denied:** inspect the port with `ls -l "$PORT"`. On
  Debian/Ubuntu it is normally owned by `dialout`; run
  `sudo usermod -aG dialout "$USER"`, then log out and back in. Follow your
  distribution's serial-port group policy elsewhere; avoid world-writable
  device permissions or running the whole toolchain as root.
- **Upload errors at high speed:** retry with `--baud 115200` for esptool, or
  `--upload-property 'upload.speed=115200'` for Arduino CLI. Recheck download
  mode and the port first.
- **Upload succeeds but no application / blank screen:** briefly reset or
  reconnect USB, then inspect serial output. Confirm the hardware is AtomS3R,
  the archive is complete, and the merged image was written at `0x0` rather
  than `0x10000`. Do not use `--force` to bypass chip/security checks.

## References

- [Project binary creation and packaging](../.github/workflows/build.yml)
- [Espressif: esptool installation](https://docs.espressif.com/projects/esptool/en/latest/esp32s3/installation.html)
- [Espressif: write/read flash and merged-image padding](https://docs.espressif.com/projects/esptool/en/latest/esp32s3/esptool/basic-commands.html)
- [Arduino CLI installation](https://docs.arduino.cc/arduino-cli/installation/)
- [Arduino CLI upload options](https://docs.arduino.cc/arduino-cli/commands-reference/arduino-cli_upload)
- [M5Stack AtomS3R hardware and download mode](https://docs.m5stack.com/en/core/AtomS3R)

[Back to sensor examples](../sensors/README.md)
