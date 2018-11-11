# ----------------------------------------------------------------
# Calculates the emittance
# Input is a pandas dataframe that has columns called x and x'.
# It returns the Emittance, sigma_w, sigma_wp and theta
# emittance [units of x * x'] [nm mrad]
# sigma_w is the semi-major axis of the ellipse [units of x]
# sigma_wp is the semi-minor axis of the ellipse [units of x']
# theta is the rotation of the ellipse [deg]
# See http://uspas.fnal.gov/materials/10MIT/Emittance.pdf
# or J. Buon, "Beam phase space and emittance".
def Calc_Emittance(df_emitt):
    sigma_x = df_emitt['x'].std(ddof=0) # \sigma_x, ddof=0 means use N as normalization
    sigma_xp = df_emitt["x'"].std(ddof=0) # \sigma_{x^\prime}
    cov_xxp = df_emitt.cov()['x']["x'"] # \sigma_x\sigma_{x^\prime}
    N = df_emitt['x'].count()
    con_xxp = cov_xxp*(N-1)/N # Use N as normalization, Pandas uses N-1
    r = df_emitt.corr(method='pearson')['x']["x'"]
    
    #emittance = sigma_x*sigma_xp*np.sqrt(1.0-r**2)
    #emittance = np.sqrt(sigma_x**2*sigma_xp**2 - cov_xxp**2)
    
    sigma_w  = np.sqrt( 0.5*(sigma_x**2 + sigma_xp**2 + np.sqrt( (sigma_x**2 - sigma_xp**2)**2 + (2.0*cov_xxp)**2 )) )
    sigma_wp = np.sqrt( 0.5*(sigma_x**2 + sigma_xp**2 - np.sqrt( (sigma_x**2 - sigma_xp**2)**2 + (2.0*cov_xxp)**2 )) )
    
    theta = 0.5*np.arctan2(2.0*cov_xxp, (sigma_x**2 - sigma_xp**2)) # in rad
    theta = theta * 180/np.pi # Convert from rad to deg
    
    #emittance = sigma_w*sigma_wp
    
    return emittance, sigma_w, sigma_wp, theta