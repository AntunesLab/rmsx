# -----------------------------------------------------------
#                VMD IDLE LOADER SCRIPT
# -----------------------------------------------------------
#
# This script waits until VMD is fully loaded (after .vmdrc)
# and then runs your main script, passing all arguments to it.

# --- CONFIGURATION ---
# ▼▼▼ EDIT THIS LINE ▼▼▼
# Set the path to your *main* script.
# IMPORTANT: Resolve relative to *this* loader script so it works from any CWD
# (e.g., when launched from Jupyter notebooks).
# You can override with environment variable RMSX_VMD_MAIN.
set script_path [file normalize [info script]]
set script_dir [file dirname $script_path]
if {[info exists env(RMSX_VMD_MAIN)] && $env(RMSX_VMD_MAIN) ne ""} {
    set main_script_path $env(RMSX_VMD_MAIN)
} else {
    set main_script_path [file join $script_dir "grid_color_scale_centered_xaxis_hotkeys.tcl"]
}
# If a relative path was supplied, resolve it against this script directory
if {[file pathtype $main_script_path] eq "relative"} {
    set main_script_path [file join $script_dir $main_script_path]
}
set main_script_path [file normalize $main_script_path]
# --- END CONFIGURATION ---


# -----------------------------------------------------------
#                VMD IDLE LOADER SCRIPT (with Delay)
# -----------------------------------------------------------
#
# This script waits until VMD is fully loaded (after .vmdrc)
# *and* a specified time has passed, then runs your main script.


# ▼▼▼ SET YOUR DELAY HERE (in milliseconds) ▼▼▼
# 1000 = 1 second
set additional_delay_ms 3000
# --- END CONFIGURATION ---


# --- FIX ---
# Set the modulation environment variable IMMEDIATELY.
# This must be done *before* VMD becomes idle, or it will be ignored.
# RMSX's main script uses:
#  - 'user' for thickness modulation
#  - 'user' for color (same values unless scaling changes)
# You can override thickness target with RMSX_VMD_THICKFIELD (user|user2|beta)
set thick_field "user"
if {[info exists env(RMSX_VMD_THICKFIELD)] && $env(RMSX_VMD_THICKFIELD) ne ""} {
    set thick_field $env(RMSX_VMD_THICKFIELD)
}
set env(VMDMODULATERIBBON) $thick_field
set env(VMDMODULATENEWTUBE) $thick_field
set env(VMDMODULATENEWCARTOON) $thick_field
puts [format {--- Set env(VMDMODULATE*) to '%s' (during init) ---} $thick_field]
# --- END FIX ---

# This procedure will run after the delay
proc run_main_script_after_delay {} {
    # Import the global command-line arguments
    global argc argv

    # Import the main script path we set above
    global main_script_path
    global env

    # Check if the file exists before trying to run it
    if {![file exists $main_script_path]} {
        puts "ERROR: Main script not found at '$main_script_path'"
        puts "Please edit wait_to_load.tcl and fix the 'main_script_path' variable."
        return
    }

    puts "--- VMD is idle and delay is over. Calling main script: $main_script_path ---"

    # Run your main script.
    source $main_script_path

    if {[info exists env(RMSX_VMD_MASK_FILE)] && $env(RMSX_VMD_MASK_FILE) ne ""} {
        if {[info procs applyMaskedTransparency] ne ""} {
            if {[catch {applyMaskedTransparency $env(RMSX_VMD_MASK_FILE)} mask_err]} {
                puts "WARNING: automatic masked transparency failed: $mask_err"
            }
        }
    }

    puts "--- Main script has finished ---"
}

# Arm the loader:
# 1) wait until VMD becomes idle (startup + .vmdrc complete)
# 2) wait an additional delay for slower systems / lots of files
puts "--- Loader armed. Main script will run after VMD is idle + $additional_delay_ms ms. ---"
after idle [list after $additional_delay_ms run_main_script_after_delay]
