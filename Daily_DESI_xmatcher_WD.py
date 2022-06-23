# After this run MWS_DA_fitter in the MWS_scripts in the DESI home directory
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import glob
import os
import sys

tmp_night_lst = sorted(glob.glob('/storage/wdplanets/DESI/Data/Daily/2*'))
last_night = int(sys.argv[1])

gaia_file = '/storage/wdplanets/phukfh/final_catalogue_v6.fits'
hdulist_gaia = fits.open(gaia_file)
gaia_wd_id = hdulist_gaia[1].data['source_id']
l_WDJ_name = hdulist_gaia[1].data['WDJ_name']
l_teff = hdulist_gaia[1].data['Teff']
l_eteff = hdulist_gaia[1].data['eTeff']
l_log_g = hdulist_gaia[1].data['log_g']
l_elog_g = hdulist_gaia[1].data['elog_g']
l_Pwd = hdulist_gaia[1].data['Pwd']
l_S_plate = hdulist_gaia[1].data['Plate']
l_S_MJD = hdulist_gaia[1].data['mjd']
l_S_fiber = hdulist_gaia[1].data['fiberID']

for k in range(len(tmp_night_lst)):
  tmp_night = tmp_night_lst[k].split('/')[-1]
  print(tmp_night)

  if int(tmp_night) >= last_night:
    

    # DAILY REDUC
    #data_dir = '/storage/astro2/phukfh/DESI/commissioning/data/{:s}/'.format(tmp_night)
    #data_dir = '/storage/astro1/phsdaj/DESI/Spectra/{:s}/'.format(tmp_night)
    #save_path = '/storage/astro2/phukfh/DESI/commissioning/SV1/'

    # EVEREST
    data_dir = '/storage/wdplanets/DESI/Data/Daily/{:s}/'.format(tmp_night)
    save_path = '/storage/astro2/phukfh/DESI/Daily/'

    meta_log = open(save_path + "csv_logs/{:s}_log.csv".format(tmp_night), "w")
    meta_log.write("#night,exposure,fiber,tileid,MJD,exptime,target_ra,target_dec,tile_ra,tile_dec,target_id,gaia_id,fiber_status,parallax,Gmag,BPmag,RPmag,H_Teff,H_Teff_e,H_log_g,H_log_g_e,WDJ_name,Pwd,SDSS_plate,SDSS_MJD,SDSS_fiber,DESI_targetbit,MWS_targetbit,StN_b,StN_r,StN_z")

    if not os.path.exists('{:s}ascii/{:s}/'.format(save_path, tmp_night)):
      os.makedirs('{:s}ascii/{:s}/'.format(save_path, tmp_night))

    if not os.path.exists('{:s}plots/{:s}/'.format(save_path, tmp_night)):
      os.makedirs('{:s}plots/{:s}/'.format(save_path, tmp_night))

    file_list = glob.glob(data_dir + "*/cframe-b?-*.fits")


    gz = 0
    if len(file_list) == 0:
      file_list = glob.glob(data_dir + "*/cframe-b?-*.fits.gz")
      gz = 1

    for i in range(len(file_list)):
      try:
        print('Openning {:s}'.format(file_list[i].split('/')[-1]))
        hdulist_DESIb = fits.open(file_list[i])
        night         = hdulist_DESIb[0].header['NIGHT']
        hdulist_DESIr = fits.open(file_list[i].replace('-b', '-r'))
        hdulist_DESIz = fits.open(file_list[i].replace('-b', '-z'))
        wave_b        = hdulist_DESIb['WAVELENGTH'].data
        wave_r        = hdulist_DESIr['WAVELENGTH'].data
        wave_z        = hdulist_DESIz['WAVELENGTH'].data
        Aflux_b       = hdulist_DESIb['FLUX'].data
        Aivar_b       = hdulist_DESIb['IVAR'].data
        Aflux_r       = hdulist_DESIr['FLUX'].data
        Aivar_r       = hdulist_DESIr['IVAR'].data
        Aflux_z       = hdulist_DESIz['FLUX'].data
        Aivar_z       = hdulist_DESIz['IVAR'].data
        tile_id       = hdulist_DESIb['FIBERMAP'].header['TILEID']
        tile_ra       = hdulist_DESIb['FIBERMAP'].header['TARGTRA']
        tile_dec      = hdulist_DESIb['FIBERMAP'].header['TARGTDEC']
        target_ra     = hdulist_DESIb['FIBERMAP'].data['TARGET_RA']
        target_dec    = hdulist_DESIb['FIBERMAP'].data['TARGET_DEC']
        target_id     = hdulist_DESIb['FIBERMAP'].data['TARGETID']
        ref_id        = hdulist_DESIb['FIBERMAP'].data['REF_ID']
        MJD           = hdulist_DESIb['FIBERMAP'].header['MJD-OBS']
        fiber_status  = hdulist_DESIb['FIBERMAP'].data['FIBERSTATUS']
        parallax      = hdulist_DESIb['FIBERMAP'].data['PARALLAX']
        exptime      = hdulist_DESIb['FIBERMAP'].data['EXPTIME']
        g_mag         = hdulist_DESIb['FIBERMAP'].data['GAIA_PHOT_G_MEAN_MAG']
        bp_mag        = hdulist_DESIb['FIBERMAP'].data['GAIA_PHOT_BP_MEAN_MAG']
        rp_mag        = hdulist_DESIb['FIBERMAP'].data['GAIA_PHOT_RP_MEAN_MAG']

        try:
          DESI_target = hdulist_DESIb['FIBERMAP'].data['DESI_TARGET']
          MWS_target  = hdulist_DESIb['FIBERMAP'].data['MWS_TARGET']
          BGS_target  = hdulist_DESIb['FIBERMAP'].data['BGS_TARGET']
        except:
          print(file_list[i], 'No SV1_XXX_TARGET in fits header')
          DESI_target = -1*np.ones(ref_id.size)
          MWS_target = -1*np.ones(ref_id.size)
          BGS_target = -1*np.ones(ref_id.size)
        # bit mask starts at 0, white dwarf is bit mask 13 (14th place)
        for j in range(ref_id.size):
          if ref_id[j] != 0:
            if ref_id[j] in gaia_wd_id:
              id = ref_id[j]
              tbit_DESI = DESI_target[j]
              tbit_MWS = MWS_target[j]
              xm = (gaia_wd_id == id)
              WDJ_name = l_WDJ_name[xm][0]
              t_id = target_id[j]
              night = str(hdulist_DESIb[0].header['NIGHT'])
              fiber = j
              big_fiber = j + int(file_list[i].split('-')[-2][-1])*500
              if gz == 0:
                exposure = file_list[i].split('-')[-1][:-5] # remove .fits
              else:
                exposure = file_list[i].split('-')[-1][:-8] # remove .fits.gz
              name = "{:s}_{:d}-{:s}-{:s}-{:04d}".format(WDJ_name,t_id,night, exposure, big_fiber)
              flux_b, ivar_b = Aflux_b[j], Aivar_b[j]
              flux_r, ivar_r = Aflux_r[j], Aivar_r[j]
              flux_z, ivar_z = Aflux_z[j], Aivar_z[j]

              StN_b = (np.sum(flux_b[(ivar_b != 0.0)] / ivar_b[(ivar_b != 0.0)]**-0.5) / flux_b[(ivar_b != 0.0)].size)
              StN_r = (np.sum(flux_r[(ivar_r != 0.0)] / ivar_r[(ivar_r != 0.0)]**-0.5) / flux_r[(ivar_r != 0.0)].size)
              StN_z = (np.sum(flux_z[(ivar_z != 0.0)] / ivar_z[(ivar_z != 0.0)]**-0.5) / flux_z[(ivar_z != 0.0)].size)

              b_sav_dat = np.transpose(np.vstack((wave_b, flux_b, ivar_b)))
              r_sav_dat = np.transpose(np.vstack((wave_r, flux_r, ivar_r)))
              z_sav_dat = np.transpose(np.vstack((wave_z, flux_z, ivar_z)))
              np.savetxt('{:s}ascii/{:s}/{:s}-b.dat'.format(save_path, night, name), b_sav_dat)
              np.savetxt('{:s}ascii/{:s}/{:s}-r.dat'.format(save_path, night, name), r_sav_dat)
              np.savetxt('{:s}ascii/{:s}/{:s}-z.dat'.format(save_path, night, name), z_sav_dat)
              plt.figure(figsize = (8,6))
              plt.plot(wave_b, flux_b, lw = 0.9)
              plt.plot(wave_r, flux_r, lw = 0.9)
              plt.plot(wave_z, flux_z, lw = 0.9)
              plt.ylim(   (    0, np.max(np.array([np.mean(flux_b) + 4*np.std(flux_b),np.mean(flux_r) + 4*np.std(flux_r),np.mean(flux_z) + 4*np.std(flux_z)]))    )    )
              plt.savefig('{:s}plots/{:s}/{:s}.pdf'.format(save_path, night, name), bbox_inches = 'tight')
              plt.close()
              line = "\n{:s},{:s},{:04d},{:d},{:.7f},{:.2f},{:.4f},{:.4f},{:.4f},{:.4f},{:d},{:d},{:d},{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:f},{:s},{:f},{:05d},{:d},{:04d},{:f},{:f},{:.2f},{:.2f},{:.2f}".format(night, exposure, big_fiber, tile_id, MJD, exptime[j], target_ra[j], target_dec[j], tile_ra, tile_dec, t_id, id, fiber_status[j], parallax[j], g_mag[j], bp_mag[j], rp_mag[j], l_teff[xm][0], l_eteff[xm][0], l_log_g[xm][0], l_elog_g[xm][0], l_WDJ_name[xm][0], l_Pwd[xm][0], l_S_plate[xm][0], l_S_MJD[xm][0], l_S_fiber[xm][0], tbit_DESI, tbit_MWS, StN_b, StN_r,StN_z)
              meta_log.write(line)
              print('DONE:', name)
              print("StN_b = {:.2f}, StN_r = {:.2f}, StN_z = {:.2f}, exptime = {:.2f}".format(StN_b, StN_r, StN_z, exptime[j]))

      except:
        print('Could not work with file {:s}'.format(file_list[i]))

