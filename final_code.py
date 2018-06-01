import numpy as np
import matplotlib.pyplot as plt

from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70., Om0=0.3)


def get_model_grid():
    """ Loads up the BPASS grid of stellar models and
    resamples it onto a coarser wavelength grid. See Step 1. """

    model_path = "data/spectra-bin-imf135_300.z020.dat"
    raw_wavelengths = np.loadtxt(model_path, usecols=0)
    raw_grid = np.loadtxt(model_path)[:,1:]

    grid = np.zeros((wavelengths.shape[0], raw_grid.shape[1]))

    for i in range(grid.shape[1]):
        grid[:,i] = np.interp(wavelengths, raw_wavelengths, raw_grid[:,i])

    grid *= (3.827*10**33)/(10**6)

    return grid


def blueshift_filters(redshift):
    """ A function that resamples filters onto the same wavelength
    basis as the model spectrum at the specified redshift. See Step 2. """

    resampled_filter_curves = []

    for filt in filter_curves:
        blueshifted_wavs = filt[:, 0]/(1 + redshift)

        resampled_filt = np.interp(wavelengths,
                                   blueshifted_wavs, filt[:, 1],
                                   left=0, right=0)

        resampled_filter_curves.append(resampled_filt)

    return resampled_filter_curves


def get_model_photometry(redshift, age_index, mass):
    """ For a row in the model grid return model
    photometry at the specified redshifts. See Step 3. """

    ssp_model = np.copy(grid[:, age_index])

    luminosity_distance = (3.086*10**24)*cosmo.luminosity_distance(redshift)

    ssp_model /= 4*np.pi*(1 + redshift)*luminosity_distance.value**2

    filter_curves_z = blueshift_filters(redshift)

    photometry = np.zeros(len(filter_curves))

    redshifted_wavs = wavelengths*(1 + redshift)

    for i in range(photometry.shape[0]):
        flux_contributions = filter_curves_z[i]*ssp_model
        photometry[i] = np.trapz(flux_contributions, x=redshifted_wavs)
        photometry[i] /= np.trapz(filter_curves_z[i], x=redshifted_wavs)

    photometry = photometry*mass

    return photometry


def load_data(row_no):
    """ Load UltraVISTA photometry from catalogue. See Step 4. """

    # load up the relevant columns from the catalogue.
    catalogue = np.loadtxt("data/UltraVISTA_catalogue.cat",
                           usecols=(0,3,4,5,6,7,8,9,10,11,12,13,14,15,
                                    16,17,18,19,20,21,22,23,24,25,26))

    # Extract the object we're interested in from the catalogue.
    fluxes = catalogue[row_no, 1:13]
    fluxerrs = catalogue[row_no, 13:25]

    # Fluxes come in erg/s/cm^2/Hz.

    # Put an error floor into the data at the 20 sigma level.
    # This is normally done to account for systematic uncertainties.
    for i in range(fluxes.shape[0]):
        if fluxerrs[i] < fluxes[i]/20:
            fluxerrs[i] = fluxes[i]/20

    # Convert from erg/s/cm^2/Hz to erg/s/cm^2/A.
    fluxes = fluxes*2.9979*10**18/eff_wavs**2
    fluxerrs = fluxerrs*2.9979*10**18/eff_wavs**2

    return fluxes, fluxerrs


# Define our basic quantities.
wavelengths = np.arange(1000., 60000., 10.)
ages = np.arange(2, 53)
ages = 10**(6+0.1*(ages-2))

grid = get_model_grid()

# Load the curves up.
filter_names = ["data/filters/CFHT_u.txt",
                "data/filters/CFHT_g.txt",
                "data/filters/CFHT_r.txt",
                "data/filters/CFHT_i+i2.txt",
                "data/filters/CFHT_z.txt",
                "data/filters/subaru_z",
                "data/filters/VISTA_Y.txt",
                "data/filters/VISTA_J.txt",
                "data/filters/VISTA_H.txt",
                "data/filters/VISTA_Ks.txt",
                "data/filters/IRAC1",
                "data/filters/IRAC2"]

filter_curves = []

for name in filter_names:
    filter_curves.append(np.loadtxt(name))

eff_wavs = np.zeros(len(filter_curves))

# Calculate the effective wavelengths of the filter curves
for i in range(len(filter_curves)):
    filt = filter_curves[i]

    wav_weights = filt[:,1]*filt[:,0]

    eff_wavs[i] = np.trapz(wav_weights, x=filt[:,0])

    eff_wavs[i] /= np.trapz(filt[:, 1], x=filt[:,0])

# Load observational data
fluxes, fluxerrs = load_data(1)

# Set grid ranges and steps
redshifts = np.arange(0.5, 1.501, 0.01)
masses = 10**np.arange(10, 11.001, 0.01)
age_indices = np.arange(51)

best_chisq = np.inf
best_point = [0.,0.,0.]

n_redshifts = redshifts.shape[0]
n_masses = masses.shape[0]
n_ages = ages.shape[0]
n_models = n_redshifts*n_masses*n_ages

# Iterate over the grid
for i in range(n_redshifts):
    for j in range(n_masses):
        for k in range(n_ages):

            model_no = k + j*n_ages + i*n_masses*n_ages

            if not model_no % 10000:
                print "Tried", model_no, "/", n_models, "models."

            redshift = redshifts[i]
            mass = masses[j]
            age_ind = age_indices[k]

            model_photometry = get_model_photometry(redshift, age_ind, mass)

            diffs = fluxes - model_photometry

            chisq = np.sum(diffs**2/fluxerrs**2)

            if chisq < best_chisq:
                best_chisq = chisq
                best_point = [i,j,k]

print "best chi-squared value:", best_chisq

print ("best parameters (z, mass, age):",
       redshifts[best_point[0]],
       np.log10(masses[best_point[1]]),
       ages[age_indices[best_point[2]]]*10**-9)

# Plot the results
best_photometry = get_model_photometry(redshifts[best_point[0]],
                                       age_indices[best_point[2]],
                                       masses[best_point[1]])

plt.figure(figsize=(15, 5))
plt.errorbar(eff_wavs, fluxes, yerr=fluxerrs, lw=1.0, linestyle=" ",
             capsize=3, capthick=1, color="black")

plt.scatter(eff_wavs, fluxes, s=75, linewidth=1,
            facecolor="blue", edgecolor="black")

plt.scatter(eff_wavs, best_photometry, s=75, facecolor="orange", zorder=10)
plt.xscale("log")
plt.xlim(1000.*(1 + redshift), 60000.*(1 + redshift))
plt.show()
