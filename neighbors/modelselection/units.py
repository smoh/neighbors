import astropy.units as u

vel_to_dist_pm = (1 * u.km/u.s).to(u.pc * u.mas/u.yr,
                                   u.dimensionless_angles()).value
dist_pm_to_vel = 1 / vel_to_dist_pm
