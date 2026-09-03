#!/usr/bin/env python3

"""
Compare DATA and MODEL_DATA from a CASA Measurement Set.

For the MS structure encountered here:

    DATA       = (4, 1, 78, 1560)
    MODEL_DATA = (4, 1, 78, 1560)
    FLAG       = (4, 1, 78, 1560)
    UVW        = (3, 78, 1560)

After removing the singleton axis:

    DATA       = (4, 78, 1560)
    MODEL_DATA = (4, 78, 1560)
    FLAG       = (4, 78, 1560)
    UVW        = (3, 78, 1560)

The script plots, averaged across the selected spectral windows:

    Top    : amplitude vs UV distance
    Bottom : phase vs UV distance

UV distance is in k-lambda.

Run from CASA 6.7.5.
"""


import numpy as np
import matplotlib.pyplot as plt

from casatools import ms, msmetadata


# ============================================================
# USER SETTINGS
# ============================================================

MS_NAME = "J1848+3219_calibrated.ms"

# Spectral windows to include.
#
# None = all spectral windows in the MS
# [0, 1, ...] = specific list of spectral windows
#
SPWS = None

# Correlation:
#
# 0 = first correlation
# 1 = second correlation
# 2 = third correlation
# 3 = fourth correlation
#
# None = average all correlations
CORR = 0

# Channel:
#
# None = all channels
# Integer = one specific channel
#
CHANNEL = None

# Bin width in k-lambda
#
# None = no binning
#
BIN_WIDTH = 50.0

# Show individual data points
SHOW_RAW_DATA = True

# Output figure
OUTPUT = "data_model_uvwave.png"


c = 299792458.0


# ============================================================
# GET SPW LIST AND CHANNEL FREQUENCIES
# ============================================================

msmd = msmetadata()

msmd.open(MS_NAME)

nspw = msmd.nspw()

msmd.done()

if SPWS is None:

    spws = list(
        range(nspw)
    )

else:

    for s in SPWS:

        if s < 0 or s >= nspw:

            raise ValueError(
                "SPW={} but MS only has {} spectral windows."
                .format(
                    s,
                    nspw
                )
            )

    spws = list(
        SPWS
    )


print()
print("SPECTRAL WINDOWS")
print("================")
print(
    "Spectral windows in MS:",
    nspw
)

print(
    "Spectral windows used:",
    spws
)


# ============================================================
# OPEN MS (ONCE)
# ============================================================

myms = ms()

myms.open(MS_NAME)


uv_all = []

data_all = []

model_all = []


for spw in spws:

    print()
    print("=" * 60)
    print("SPW", spw)
    print("=" * 60)

    # --------------------------------------------------------
    # Channel frequencies for this SPW
    # --------------------------------------------------------

    msmd = msmetadata()

    msmd.open(MS_NAME)

    freqs = np.asarray(
        msmd.chanfreqs(spw)
    )

    msmd.done()


    print()
    print("Channel frequencies:")
    print(
        "  Number of channels in SPW header:",
        len(freqs)
    )

    print(
        "  Frequency range:",
        freqs.min() / 1e9,
        "-",
        freqs.max() / 1e9,
        "GHz"
    )

    # --------------------------------------------------------
    # Select SPW and get data
    # --------------------------------------------------------

    myms.reset()

    myms.msselect({
        "spw": str(spw)
    })

    dat = myms.getdata(
        items=[
            "corrected_data",
            "model_data",
            "flag",
            "uvw"
        ],
        ifraxis=True
    )

    # --------------------------------------------------------
    # Guards: empty selection for this SPW
    # --------------------------------------------------------

    if "corrected_data" not in dat:

        print()
        print(
            "WARNING: no data for SPW",
            spw,
            "-- skipping."
        )

        continue

    # --------------------------------------------------------
    # Get arrays
    # --------------------------------------------------------

    DATA = np.squeeze(
        np.asarray(
            dat["corrected_data"]
        )
    )

    MODEL = np.squeeze(
        np.asarray(
            dat["model_data"]
        )
    )

    FLAG = np.squeeze(
        np.asarray(
            dat["flag"]
        )
    )

    UVW = np.asarray(
        dat["uvw"]
    )


    print()
    print("DATA       :", DATA.shape)
    print("MODEL_DATA :", MODEL.shape)
    print("FLAG       :", FLAG.shape)
    print("UVW        :", UVW.shape)


    # --------------------------------------------------------
    # Check expected dimensions
    # --------------------------------------------------------

    if DATA.ndim != 3:

        raise RuntimeError(
            "Unexpected DATA shape after squeeze: {}".format(
                DATA.shape
            )
        )


    if UVW.ndim != 3:

        raise RuntimeError(
            "Unexpected UVW shape: {}".format(
                UVW.shape
            )
        )


    npol = DATA.shape[0]

    nchan = DATA.shape[1]

    nrow = DATA.shape[2]


    if UVW.shape[0] != 3:

        raise RuntimeError(
            "Expected UVW first dimension to be 3. "
            "Got {}".format(UVW.shape)
        )


    if UVW.shape[1] != nchan:

        raise RuntimeError(
            "UVW channel dimension does not match DATA."
            "\nUVW  = {}".format(UVW.shape)
            + "\nDATA = {}".format(DATA.shape)
        )


    if UVW.shape[2] != nrow:

        raise RuntimeError(
            "UVW row dimension does not match DATA."
            "\nUVW  = {}".format(UVW.shape)
            + "\nDATA = {}".format(DATA.shape)
        )


    # --------------------------------------------------------
    # Select correlation
    # --------------------------------------------------------

    if CORR is None:

        print()
        print("Averaging all correlations.")

        DATA = np.mean(
            DATA,
            axis=0
        )

        MODEL = np.mean(
            MODEL,
            axis=0
        )

        FLAG = np.any(
            FLAG,
            axis=0
        )

    else:

        if CORR >= npol:

            raise ValueError(
                "CORR={} but DATA only has {} correlations."
                .format(
                    CORR,
                    npol
                )
            )

        print()
        print(
            "Using correlation:",
            CORR
        )

        DATA = DATA[
            CORR,
            :,
            :
        ]

        MODEL = MODEL[
            CORR,
            :,
            :
        ]

        FLAG = FLAG[
            CORR,
            :,
            :
        ]


    # --------------------------------------------------------
    # Channel selection
    # --------------------------------------------------------

    if CHANNEL is None:

        channels = np.arange(
            nchan
        )

    else:

        if CHANNEL >= nchan:

            raise ValueError(
                "CHANNEL={} but only {} channels."
                .format(
                    CHANNEL,
                    nchan
                )
            )

        channels = np.array(
            [CHANNEL]
        )


    # --------------------------------------------------------
    # UV distance
    # --------------------------------------------------------

    #
    # IMPORTANT:
    #
    # Your UVW has shape:
    #
    #       (3, channel, row)
    #
    # Therefore:
    #
    #       U = UVW[0,:,:]
    #       V = UVW[1,:,:]
    #
    # gives:
    #
    #       (channel, row)
    #
    # which matches DATA.
    #

    U = UVW[0, :, :]

    V = UVW[1, :, :]


    # UV distance in metres

    uv_m = np.sqrt(
        U**2 + V**2
    )


    # wavelengths:
    #
    # shape = (channel,)

    wavelengths = (
        c / freqs[:nchan]
    )


    # Expand wavelength to:
    #
    #     (channel, 1)
    #
    # so it broadcasts correctly against
    # uv_m = (channel, row)

    uv_lambda = (
        uv_m
        /
        wavelengths[:, None]
    )


    # Convert to k-lambda

    uv_klambda = (
        uv_lambda / 1000.0
    )


    # --------------------------------------------------------
    # Collect valid data for this SPW
    # --------------------------------------------------------

    for ch in channels:

        uv_ch = uv_klambda[ch, :]

        data_ch = DATA[ch, :]

        model_ch = MODEL[ch, :]

        flag_ch = FLAG[ch, :]


        valid = (
            (~flag_ch)
            &
            np.isfinite(uv_ch)
            &
            np.isfinite(data_ch.real)
            &
            np.isfinite(data_ch.imag)
            &
            np.isfinite(model_ch.real)
            &
            np.isfinite(model_ch.imag)
        )


        uv_all.append(
            uv_ch[valid]
        )

        data_all.append(
            data_ch[valid]
        )

        model_all.append(
            model_ch[valid]
        )


myms.close()


# ============================================================
# COMBINE ALL SPWS / CHANNELS
# ============================================================

uv = np.concatenate(
    uv_all
)

data_vis = np.concatenate(
    data_all
)

model_vis = np.concatenate(
    model_all
)


print()
print("VALID DATA")
print("==========")
print(
    "Number of visibilities:",
    len(uv)
)


# ============================================================
# AMPLITUDE
# ============================================================

data_amp = np.abs(
    data_vis
)

model_amp = np.abs(
    model_vis
)


# ============================================================
# PHASE
# ============================================================

data_phase = np.angle(
    data_vis,
    deg=True
)

model_phase = np.angle(
    model_vis,
    deg=True
)


# ============================================================
# BINNING FUNCTION
# ============================================================

def bin_median(
    x,
    y,
    width
):

    if width is None:

        return (
            x,
            y
        )


    xmin = np.nanmin(
        x
    )

    xmax = np.nanmax(
        x
    )


    bins = np.arange(
        xmin,
        xmax + width,
        width
    )


    xbin = []

    ybin = []


    for i in range(
        len(bins) - 1
    ):

        mask = (
            (x >= bins[i])
            &
            (x < bins[i + 1])
        )


        if not np.any(mask):

            continue


        xbin.append(
            0.5 * (
                bins[i]
                +
                bins[i + 1]
            )
        )


        ybin.append(
            np.nanmedian(
                y[mask]
            )
        )


    return (
        np.asarray(xbin),
        np.asarray(ybin)
    )


# ============================================================
# BIN DATA
# ============================================================

uv_d_amp, amp_d = bin_median(
    uv,
    data_amp,
    BIN_WIDTH
)

uv_m_amp, amp_m = bin_median(
    uv,
    model_amp,
    BIN_WIDTH
)


uv_d_phase, phase_d = bin_median(
    uv,
    data_phase,
    BIN_WIDTH
)

uv_m_phase, phase_m = bin_median(
    uv,
    model_phase,
    BIN_WIDTH
)


# ============================================================
# MAKE FIGURE
# ============================================================

fig, axes = plt.subplots(
    2,
    1,
    figsize=(10, 9),
    sharex=True
)


# ============================================================
# AMPLITUDE PANEL
# ============================================================

ax = axes[0]

if SHOW_RAW_DATA:

    ax.scatter(
        uv,
        data_amp,
        s=3,
        alpha=0.15,
        label="DATA"
    )

# Binned DATA
ax.plot(
    uv_d_amp,
    amp_d,
    ".",
    markersize=5,
    label="DATA median"
)

# MODEL_DATA -- points only, NO connecting line
ax.plot(
    uv_m_amp,
    amp_m,
    ".",
    markersize=7,
    label="MODEL_DATA"
)

ax.set_ylabel("Amplitude")
ax.set_title("Visibility amplitude vs UV distance")

ax.legend()
ax.grid(alpha=0.25)


# ============================================================
# PHASE PANEL
# ============================================================

ax = axes[1]

if SHOW_RAW_DATA:

    ax.scatter(
        uv,
        data_phase,
        s=3,
        alpha=0.15,
        label="DATA"
    )

# Binned DATA
ax.plot(
    uv_d_phase,
    phase_d,
    ".",
    markersize=5,
    label="DATA median"
)

# MODEL_DATA -- points only, NO connecting line
ax.plot(
    uv_m_phase,
    phase_m,
    ".",
    markersize=7,
    label="MODEL_DATA"
)

ax.set_xlabel(
    r"UV distance (k$\lambda$)"
)

ax.set_ylabel("Phase (deg)")

ax.set_ylim(-180, 180)

ax.legend()
ax.grid(alpha=0.25)


# ============================================================
# SAVE
# ============================================================

plt.tight_layout()


plt.savefig(
    OUTPUT,
    dpi=200,
    bbox_inches="tight"
)


plt.show()


print()
print("Figure saved as:")
print(OUTPUT)
