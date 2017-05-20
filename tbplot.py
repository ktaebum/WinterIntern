import matplotlib.pyplot as plt
import numpy as np


def plot_histo(disk, num):
    print "\n\nMaking Histograms..."
    fig, ax = plt.subplots(1, 3, figsize=(16, 9))
    for i in range(3):
        if i == 0 or i == 1:
            ax[i].hist(disk[:, i], bins=100, facecolor='white', range=(-30, 30))
        else:
            ax[i].hist(disk[:, i], bins=100, facecolor='white', range=(-10, 10))
    print "Finished Making Histogram"
    print "Saving Histogram image..."
    fig.tight_layout()
    fig.savefig('histo_{:03d}.png'.format(num))
    print 'Finished Saving'
    plt.gcf().clear()
    plt.close()


def plot_particle(disk, num):
    print '\n\nPlotting Particles Image...'
    plt.figure(figsize=(12, 12))
    plt.plot(disk[0], disk[1], '*', color='black')
    plt.xlim(-30, 30)
    plt.ylim(-30, 30)
    plt.title('X-Y Particles')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('position_{:03d}.png'.format(num))
    plt.gcf().clear()
    plt.close()
    print 'Finished Plotting and Saving'


def plot_contour(grid, rho, rho_bar, num, name):
    print "\n\nNow Plotting Contour Graph With Surface Density..."
    plt.figure(figsize=(15, 12))
    plt.set_cmap('magma')
    plt.rc('font', size=22)
    # plt.contourf(grid, grid, rho, levels=np.mgrid[np.power(np.e, -13):rho_bar:256j])
    plt.contourf(grid, grid, np.log(rho), levels=np.mgrid[-11:np.log(rho_bar):256j])
    plt.title('Density (log)', fontsize=30)
    plt.xlim(-8, 8)
    plt.xticks(np.mgrid[-8:8:5j])
    plt.ylim(-8, 8)
    plt.yticks(np.mgrid[-8:8:5j])
    plt.colorbar()
    # plt.xlabel('x [kpc]', fontsize=30)
    # plt.ylabel('z [kpc]', fontsize=30)
    plt.figtext(0.15, 0.83, 'Snapshot = {:d}'.format(num), fontsize=28, color='black')
    # plt.figtext(0.15, 0.83, 'T = {:4.3f}Gyr'.format(num / 1000.), fontsize=28, color='black')
    plt.savefig(str(name) + '.png')
    # plt.savefig('new_surface_density_{:03d}.png'.format(num))
    plt.gcf().clear()
    plt.close()
    print "Finished Plotting Contour Graph\n"


def plot_bar_fit(grid, x_sli, y_sli, z_sli, rho_bar, a, b, vector, num):
    ylim = 0.75
    ygrid = 5
    xlim = 5
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    ax[0][0].plot(grid, x_sli)
    ax[0][0].plot(a, rho_bar * np.power(1 - np.power(a / vector[1], 2), vector[0]))
    ax[0][0].set_xlim(-xlim, xlim)
    ax[0][0].set_xticks([])
    ax[0][0].set_ylim(0, ylim)
    ax[0][0].set_yticks(np.mgrid[0:ylim:complex(0, ygrid)])
    ax[0][0].tick_params(axis='x', labelsize=15)
    ax[0][0].tick_params(axis='y', labelsize=15)
    ax[0][0].set_ylabel(r'$\rho$', fontsize=20)
    ax[0][0].text(-4.8, ylim * 0.73, 'Snapshot = {:d}\nn = {:5.4f}\nmajor = {:5.4f}\nminor = {:5.4f}'
                  .format(num, vector[0], vector[1], vector[2]), fontsize=15)
    ax[0][0].set_title('Bar Fitting', fontsize=25, loc='left')
    ax[0][1].plot(grid, y_sli)
    ax[0][1].plot(b, rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0]))
    ax[0][1].set_xlim(-xlim, xlim)
    ax[0][1].set_xticks([])
    ax[0][1].set_ylim(0, ylim)
    ax[0][1].set_yticks([])
    ax[0][1].tick_params(axis='x', labelsize=15)
    ax[0][1].tick_params(axis='y', labelsize=15)
    ax[0][2].plot(grid, z_sli, label='Real')
    ax[0][2].plot(b, rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0]), label='Fit')
    ax[0][2].set_xlim(-xlim, xlim)
    ax[0][2].set_xticks([])
    ax[0][2].set_ylim(0, ylim)
    ax[0][2].set_yticks([])
    ax[0][2].tick_params(axis='x', labelsize=15)
    ax[0][2].tick_params(axis='y', labelsize=15)
    ax[0][2].legend()
    xlim_u = 5
    xlim_l = -5
    ylim_u = np.log(ylim)
    ylim_l = -5
    x_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    y_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    z_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    x_sli = np.log(x_sli)
    y_sli = np.log(y_sli)
    z_sli = np.log(z_sli)
    ax[1][0].plot(grid, x_sli)
    ax[1][0].plot(a, np.log(rho_bar * np.power(1 - np.power(a / vector[1], 2), vector[0])))
    ax[1][0].set_xlim(xlim_l, xlim_u)
    ax[1][0].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[1][0].set_ylim(ylim_l, ylim_u)
    ax[1][0].set_yticks(np.mgrid[ylim_l:ylim_u:complex(0, ygrid)])
    ax[1][0].tick_params(axis='x', labelsize=15)
    ax[1][0].tick_params(axis='y', labelsize=15)
    ax[1][0].set_xlabel('x[kpc]', fontsize=20)
    ax[1][0].set_ylabel(r'$log(\rho)$', fontsize=20)
    ax[1][1].plot(grid, y_sli)
    ax[1][1].plot(b, np.log(rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0])))
    ax[1][1].set_xlim(xlim_l, xlim_u)
    ax[1][1].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[1][1].set_ylim(ylim_l, ylim_u)
    ax[1][1].tick_params(axis='x', labelsize=15)
    ax[1][1].set_yticks([])
    ax[1][1].tick_params(axis='y', labelsize=15)
    ax[1][1].set_xlabel('y[kpc]', fontsize=20)
    ax[1][2].plot(grid, z_sli)
    ax[1][2].plot(b, np.log(rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0])))
    ax[1][2].set_xlim(xlim_l, xlim_u)
    ax[1][2].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[1][2].set_ylim(ylim_l, ylim_u)
    ax[1][2].set_yticks([])
    ax[1][2].tick_params(axis='x', labelsize=15)
    ax[1][2].tick_params(axis='y', labelsize=15)
    ax[1][2].set_xlabel('z[kpc]', fontsize=20)
    fig.tight_layout()
    plt.savefig('compare_{:03d}.png'.format(num))
    plt.close()


def plot_bar_fit_log(grid, x_sli, y_sli, z_sli, rho_bar, a, b, vector, num):
    xlim_u = 5
    xlim_l = -5
    ylim_u = np.log(np.e)
    ylim_l = -5
    x_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    y_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    z_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    x_sli = np.log(x_sli)
    y_sli = np.log(y_sli)
    z_sli = np.log(z_sli)
    ygrid = 5
    fig, ax = plt.subplots(1, 3, figsize=(16, 9))
    ax[0].plot(grid, x_sli)
    ax[0].plot(a, np.log(rho_bar * np.power(1 - np.power(a / vector[1], 2), vector[0])))
    ax[0].set_xlim(xlim_l, xlim_u)
    ax[0].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[0].set_ylim(ylim_l, ylim_u)
    ax[0].set_yticks(np.mgrid[ylim_l:ylim_u:complex(0, ygrid)])
    ax[0].tick_params(axis='x', labelsize=15)
    ax[0].tick_params(axis='y', labelsize=15)
    ax[0].set_xlabel('x[kpc]', fontsize=20)
    ax[0].set_ylabel(r'$log(\rho)$', fontsize=20)
    ax[0].text(-4.8, 0., 'Snapshot = {:d}\nn = {:5.4f}\nmajor = {:5.4f}\nminor = {:5.4f}'
               .format(num, vector[0], vector[1], vector[2]), fontsize=20)
    ax[1].plot(grid, y_sli)
    ax[1].plot(b, np.log(rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0])))
    ax[1].set_xlim(xlim_l, xlim_u)
    ax[1].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[1].set_ylim(ylim_l, ylim_u)
    ax[1].set_yticks([])
    ax[1].tick_params(axis='x', labelsize=15)
    ax[1].tick_params(axis='y', labelsize=15)
    ax[1].set_xlabel('y[kpc]', fontsize=20)
    ax[2].plot(grid, z_sli)
    ax[2].plot(b, np.log(rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0])))
    ax[2].set_xlim(xlim_l, xlim_u)
    ax[2].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[2].set_ylim(ylim_l, ylim_u)
    ax[2].set_yticks([])
    ax[2].tick_params(axis='x', labelsize=15)
    ax[2].tick_params(axis='y', labelsize=15)
    ax[2].set_xlabel('z[kpc]', fontsize=20)
    fig.tight_layout()
    plt.savefig('compare_{:03d}_log.png'.format(num))


def plot_fitcompare(a, candidate, theta, major, minor, chi, num):
    slc = np.where(np.power(a[1] * candidate[0, :], 2) - 4 * a[2] * (a[3] + a[0] * np.power(candidate[0, :], 2)) >= 0)
    y_f1, y_f2 = fit_ellipse(candidate[0, :], a, slc[0])
    text = '$\phi = {:.4f}rad, {:.4f}^\degree$\n'.format(theta, theta * 180. / np.pi)
    text += '$Major = {:.4f}$\n$Minor = {:.4f}$\n'.format(max(major, minor), min(major, minor))
    text += r'$D(\vec u)/N = {:.4e}$'.format(chi)
    print '\n\nPlotting FitCompare Image...'
    plt.figure(figsize=(12, 12))
    plt.rc('font', size=22)
    plt.plot(candidate[0, :], candidate[1, :], '^', color='red', label='RealData')
    plt.plot(candidate[0, slc[0]], y_f1, '.', color='blue', label='FittedData')
    plt.plot(candidate[0, slc[0]], y_f2, '.', color='blue')
    plt.title('Compare Fitting', fontsize=30)
    plt.xlim(-6, 6)
    plt.ylim(-6, 6)
    plt.xlabel('X', fontsize=30)
    plt.ylabel('Y', fontsize=30)
    plt.legend(fontsize=30)
    plt.figtext(0.15, 0.14, text, fontsize=28, fontstyle='italic')
    plt.figtext(0.15, 0.82, 'T = {:4.3f}Gyr'.format(num / 1000.), fontsize=28)
    plt.savefig('fit_{:03d}.png'.format(num))
    plt.gcf().clear()
    plt.close()
    print 'Finished Plotting'


def fit_ellipse(x, u, slc):
    insqrt = np.power(u[1] * x[slc], 2) - 4 * u[2] * (u[3] + u[0] * np.power(x[slc], 2))
    yfit_p = (-u[1] * x[slc] + np.sqrt(insqrt)) / (2. * u[2])
    yfit_m = (-u[1] * x[slc] - np.sqrt(insqrt)) / (2. * u[2])
    return yfit_p, yfit_m


def plot_bulge_fit(grid, rho, rho_bar, radius, num):
    fig, ax = plt.subplots(1, 3, sharey='row', figsize=(16, 9))
    xlim_u = 5
    xlim_l = -5
    ylim_u = 1
    ylim_l = 0
    ygrid = 5
    n = len(grid)
    x_sli = rho[n / 2, n / 2, :]
    y_sli = rho[n / 2, :, n / 2]
    z_sli = rho[:, n / 2, n / 2]
    x = np.linspace(-radius, radius, 100)

    real, = ax[0].plot(grid, x_sli, label='Real')
    fit, = ax[0].plot(x, rho_bar * np.power(1 + np.power(x / radius, 2), -1.5), label='Fit')
    ax[0].set_xlim(xlim_l, xlim_u)
    ax[0].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[0].set_ylim(ylim_l, ylim_u)
    ax[0].set_yticks(np.mgrid[ylim_l:ylim_u:complex(0, ygrid)])
    ax[0].tick_params(axis='x', labelsize=15)
    ax[0].tick_params(axis='y', labelsize=15)
    ax[0].set_xlabel('x[kpc]', fontsize=20)
    ax[0].set_ylabel(r'$\rho$', fontsize=20)
    ax[0].set_title('Bulge Fit', loc='left')
    ax[0].text(-4.8, ylim_u * 0.92, 'Snapshot = {:d}\nradius = {:5.4f}'.format(num, radius), fontsize=20)
    ax[1].plot(grid, y_sli)
    ax[1].plot(x, rho_bar * np.power(1 + np.power(x / radius, 2), -1.5))
    ax[1].set_xlim(xlim_l, xlim_u)
    ax[1].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[1].set_ylim(ylim_l, ylim_u)
    ax[1].set_yticks(np.mgrid[ylim_l:ylim_u:complex(0, ygrid)])
    ax[1].tick_params(axis='x', labelsize=15)
    ax[1].tick_params(axis='y', labelsize=15)
    ax[1].set_xlabel('y[kpc]', fontsize=20)

    ax[2].plot(grid, z_sli)
    ax[2].plot(x, rho_bar * np.power(1 + np.power(x / radius, 2), -1.5))
    ax[2].set_xlim(xlim_l, xlim_u)
    ax[2].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[2].set_ylim(ylim_l, ylim_u)
    ax[2].set_yticks(np.mgrid[ylim_l:ylim_u:complex(0, ygrid)])
    ax[2].tick_params(axis='x', labelsize=15)
    ax[2].tick_params(axis='y', labelsize=15)
    ax[2].set_xlabel('z[kpc]', fontsize=20)
    plt.legend(handles=[real, fit], loc='upper right')
    fig.tight_layout()
    plt.savefig('bulgefit_{:03d}.png'.format(num))
    plt.close()


def plot_cont_fit(grid, xy_rho, rho_bar, candidate, a, theta, major, minor, chi_fit, num):
    fig, ax = plt.subplots(1, 2, figsize=(16, 9))
    plt.set_cmap('magma')
    ax[0].contourf(grid, grid, xy_rho, levels=np.mgrid[np.power(np.e, -13):rho_bar:256j], fontsize=22)
    # ax[0].contourf(grid, grid, np.log(xy_rho), levels=np.mgrid[-13:np.log(rho_bar):256j], fontsize=22)
    ax[0].set_title('Density (linear)', fontsize=30)
    ax[0].set_xlim(-8, 8)
    ax[0].set_ylim(-8, 8)
    ax[0].set_xticks(np.mgrid[-8:8:5j])
    ax[0].tick_params(axis='x', labelsize=20)
    ax[0].tick_params(axis='y', labelsize=20)
    ax[0].set_yticks(np.mgrid[-8:8:5j])
    ax[0].set_xlabel('X[kpc]', fontsize=25)
    ax[0].set_ylabel('Y[kpc]', fontsize=25)
    ax[0].text(-7.5, 7, 'Snapshot = {:d}'.format(num), fontsize=25, color='white')
    slc = np.where(
        np.power(a[1] * candidate[0, :], 2) - 4 * a[2] * (a[3] + a[0] * np.power(candidate[0, :], 2)) >= 0)
    y_f1, y_f2 = fit_ellipse(candidate[0, :], a, slc[0])
    text = '$\phi = {:.4f}rad, {:.4f}^\degree$\n'.format(theta, theta * 180. / np.pi)
    text += '$Major = {:.4f}$\n$Minor = {:.4f}$\n'.format(max(major, minor), min(major, minor))
    text += r'$D(\vec u)/N = {:.4e}$'.format(chi_fit)
    ax[1].plot(candidate[0, :], candidate[1, :], '^', color='red', label='RealData')
    ax[1].plot(candidate[0, slc[0]], y_f1, '.', color='blue', label='FittedData')
    ax[1].plot(candidate[0, slc[0]], y_f2, '.', color='blue')
    ax[1].set_title('Ellipse Fitting', fontsize=30)
    ax[1].set_xlim(-8, 8)
    ax[1].set_ylim(-8, 8)
    ax[1].set_xticks(np.mgrid[-8:8:5j])
    ax[1].tick_params(axis='x', labelsize=20)
    ax[1].tick_params(axis='y', labelsize=20)
    # ax[1].set_yticks(np.mgrid[-8:8:9j])
    ax[1].set_yticks([])
    ax[1].set_xlabel('X[kpc]', fontsize=25)
    # ax[1].set_ylabel('Y[kpc]', fontsize=25)
    ax[1].legend(fontsize=25)
    ax[1].text(-7.5, -7.5, text, fontsize=25, fontstyle='italic')
    # ax[1].text(0.15, 0.82, 'T = {:4.3f}Gyr'.format(num / 1000.), fontsize=28)
    fig.tight_layout()
    plt.savefig('{:d}'.format(num))
    plt.gcf().clear()
    plt.close()


def plot_bar_fit22(grid, x_sli, y_sli, z_sli, rho_bar, a, b, vector, num):
    ylim = 0.7
    ygrid = 5
    xlim = 5
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    ax[0][0].plot(grid, x_sli)
    ax[0][0].plot(a, rho_bar * np.power(1 - np.power(a / vector[1], 2), vector[0]))
    ax[0][0].set_xlim(-xlim, xlim)
    ax[0][0].set_xticks([])
    ax[0][0].set_ylim(0, ylim)
    ax[0][0].set_yticks(np.mgrid[0:ylim:complex(0, ygrid)])
    ax[0][0].tick_params(axis='x', labelsize=15)
    ax[0][0].tick_params(axis='y', labelsize=15)
    ax[0][0].set_ylabel(r'$\rho$', fontsize=20)
    ax[0][0].text(-4.8, ylim * 0.76, 'Snapshot = {:d}\nn = {:5.4f}\nmajor = {:5.4f}\nminor = {:5.4f}'
                  .format(num, vector[0], vector[1], vector[2]), fontsize=15)
    ax[0][1].plot(grid, y_sli)
    ax[0][1].plot(b, rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0]))
    ax[0][1].set_xlim(-xlim, xlim)
    ax[0][1].set_xticks([])
    ax[0][1].set_ylim(0, ylim)
    ax[0][1].set_yticks([])
    ax[0][1].tick_params(axis='x', labelsize=15)
    ax[0][1].tick_params(axis='y', labelsize=15)
    ax[0][2].plot(grid, z_sli, label='Real')
    ax[0][2].plot(b, rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0]), label='Fit')
    ax[0][2].set_xlim(-xlim, xlim)
    ax[0][2].set_xticks([])
    ax[0][2].set_ylim(0, ylim)
    ax[0][2].set_yticks([])
    ax[0][2].tick_params(axis='x', labelsize=15)
    ax[0][2].tick_params(axis='y', labelsize=15)
    xlim_u = 5
    xlim_l = -5
    ylim_u = np.log(np.e)
    ylim_l = -5
    x_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    y_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    z_sli[x_sli < np.power(np.e, -11)] = np.power(np.e, -11)
    x_sli = np.log(x_sli)
    y_sli = np.log(y_sli)
    z_sli = np.log(z_sli)
    ax[1][0].plot(grid, x_sli)
    ax[1][0].plot(a, np.log(rho_bar * np.power(1 - np.power(a / vector[1], 2), vector[0])))
    ax[1][0].set_xlim(xlim_l, xlim_u)
    ax[1][0].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[1][0].set_ylim(ylim_l, ylim_u)
    ax[1][0].set_yticks(np.mgrid[ylim_l:ylim_u:complex(0, ygrid)])
    ax[1][0].tick_params(axis='x', labelsize=15)
    ax[1][0].tick_params(axis='y', labelsize=15)
    ax[1][0].set_xlabel('x[kpc]', fontsize=20)
    ax[1][0].set_ylabel(r'$log(\rho)$', fontsize=20)
    ax[1][1].plot(grid, y_sli)
    ax[1][1].plot(b, np.log(rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0])))
    ax[1][1].set_xlim(xlim_l, xlim_u)
    ax[1][1].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[1][1].set_ylim(ylim_l, ylim_u)
    ax[1][1].set_yticks([])
    ax[1][1].tick_params(axis='x', labelsize=15)
    ax[1][1].tick_params(axis='y', labelsize=15)
    ax[1][1].set_xlabel('y[kpc]', fontsize=20)
    ax[1][2].plot(grid, z_sli)
    ax[1][2].plot(b, np.log(rho_bar * np.power(1 - np.power(b / vector[2], 2), vector[0])))
    ax[1][2].set_xlim(xlim_l, xlim_u)
    ax[1][2].set_xticks(np.array([xlim_l, 0, xlim_u]))
    ax[1][2].set_ylim(ylim_l, ylim_u)
    ax[1][2].set_yticks([])
    ax[1][2].tick_params(axis='x', labelsize=15)
    ax[1][2].tick_params(axis='y', labelsize=15)
    ax[1][2].set_xlabel('z[kpc]', fontsize=20)
    fig.tight_layout()
    plt.savefig('compare_{:03d}.png'.format(num))


def plot_zvel(grid, vel, num):
    vel = abs(vel)
    mean = np.mean(vel)
    plt.figure(figsize=(12, 12))
    plt.plot(grid, vel)
    plt.savefig('v_vel_{:d}'.format(num))
    plt.close()
