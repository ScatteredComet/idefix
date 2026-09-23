import marimo

__generated_with = "0.24.2"
app = marimo.App()


@app.cell
def _():
    import os, sys

    sys.path.append(
        os.path.join(os.environ.get("IDEFIX_DIR", "../../../../"), "pytools")
    )

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import dump_io as io

    from astropy import constants as const, units as u

    # '%matplotlib inline' command supported automatically in marimo
    return LogNorm, const, io, np, plt, u


@app.cell
def _(PrimitiveCenter, io, np):
    output_directory = "out/"

    # vc is the "main cell-centered primitive variables index"
    # rho is the ???

    def load(dirc=output_directory, time_slice_index=0, var=PrimitiveCenter.Density, slicer=np.s_[:]):
        return (
            io.readDump(f"{output_directory}/dump.{time_slice_index:04d}.dmp")
            .data[var]
            .transpose(2, 1, 0)[slicer][0]
        )

    def load3D(dirc=output_directory, time_slice_index=0, var=PrimitiveCenter.Density, slicer=np.s_[:]):
        return (
            io.readDump(f"{output_directory}/dump.{time_slice_index:04d}.dmp")
            .data[var]
            .transpose(2, 1, 0)[slicer]
        )

    return load, load3D, output_directory


@app.cell
def _(const, io, np, output_directory, u):
    dmp = io.readDump(f"{output_directory}/dump.0000.dmp")
    r, th, _phi = dmp.x1, dmp.x2, dmp.x3

    # define (reference) units of length, mass, velocity, temperature and density
    # these are astropy units and will be converted later (e.g. uT to Kelvin) to give conversion factors
    # i assume these are defined in the cpp code somewhere, idk if in the setup.cpp or in idefix?
    uL = const.au.cgs
    uM = const.M_sun.cgs
    uv = np.sqrt(const.G.cgs * uM / uL)
    uT = (
        (u.g.cgs / u.mol.cgs) * uv**2 / (const.k_B.cgs * const.N_A.cgs)
    ).decompose()
    ud = uM / uL**3

    mu = 2.353 # standard mostly hydrogen + some helium assumption
    return dmp, mu, r, th, uT, ud


@app.cell
def _():
    # define strings used within dmp.data.keys()

    # V or W are used for primitive quantities (e.g. density)
    # U is used for conserved quantities (e.g. density per unit volume, momentum)
    # c means evaluated at the center of a cell

    from enum import StrEnum

    class PrimitiveCenter(StrEnum):
        """
        Primitive aka NOT per unit volume (?)
        """
        Pressure = "Vc-PRS"
        Density = "Vc-RHO"

    class Misc(StrEnum):
        RadiationEnergyDensity = "ERAD"


    return Misc, PrimitiveCenter


@app.cell
def _(LogNorm, Misc, PrimitiveCenter, load, mu, np, plt, r, th, uT, ud):
    _fig, _ax = plt.subplots(ncols=2, dpi=150, figsize=(8, 3))
    snapshot = 10
    rho = load(time_slice_index=snapshot, var=PrimitiveCenter.Density)
    P = load(time_slice_index=snapshot, var=PrimitiveCenter.Pressure)
    Er = load(time_slice_index=snapshot, var=Misc.RadiationEnergyDensity)
    _T = mu * (P / rho)

    # NOTE: this is a 2D array whose second index indexes radius
    z = r * np.cos(th[:, None]) # None converts th[] to a 2D array (?)

    im = _ax[0].pcolormesh(
        r,
        z,
        rho * ud.to_value("g/cm**3"),
        norm=LogNorm(vmin=1e-10 * rho.max() * ud.to_value("g/cm**3")),
    )
    _ax[0].set_xlabel("Radius [AU]")
    _ax[0].set_ylabel("Z/R")
    _fig.colorbar(im, ax=_ax[0], label="Density [cgs]")
    im = _ax[1].pcolormesh(r, z, _T * uT.to_value("K"), norm=LogNorm())
    _ax[1].set_xlabel("Radius [AU]")
    _ax[1].set_ylabel("Z/R")
    _fig.colorbar(im, ax=_ax[1], label="Temperature [K]")
    _fig.tight_layout()
    plt.show()
    return snapshot, z


@app.cell
def _(PrimitiveCenter, load3D, mu, plt, r, snapshot, uT):
    _fig, _ax = plt.subplots(ncols=1, dpi=150)
    _T = mu * (load3D(time_slice_index=snapshot, var=PrimitiveCenter.Pressure) / load3D(time_slice_index=snapshot, var=PrimitiveCenter.Density))
    for _i in range(_T.shape[0]):
        _ax.loglog(r, _T[_i, -1, :] * uT.to_value("K"), ls="-")
        _ax.loglog(r, _T[_i, 0, :] * uT.to_value("K"), ls="--")
    _ax.loglog(r, 350 * r ** (-0.5), c="k", ls="--")
    _ax.set_xlabel("Radius [AU]")
    _ax.set_ylabel("T [K]")
    _fig.tight_layout()
    plt.show()
    return


@app.cell
def _(PrimitiveCenter, load3D, mu, plt, snapshot, uT, z):
    _fig, _ax = plt.subplots(ncols=1, dpi=150)
    _T = mu * (load3D(time_slice_index=snapshot, var=PrimitiveCenter.Pressure) / load3D(time_slice_index=snapshot, var=PrimitiveCenter.Density))
    for _i in [10, 30, 60]: # radii indices
        # z is a 2D array whose 2nd index indexes radius.
        # its first index indexes z
        # so _T is a 3D array whose 3rd index indexes radius, 2nd index indexes z, and 1st index indexes ???
        _ax.semilogy(z[:, _i], _T[0, :, _i] * uT.to_value("K"), ls="-")
    _ax.set_xlabel("Z [AU]")
    _ax.set_ylabel("T [K]")
    _fig.tight_layout()
    plt.show()
    return


@app.cell
def _(dmp):
    print(dmp.data.keys())
    return


if __name__ == "__main__":
    app.run()
