import csv
import os
import sys
from tbcal import *
from tbfit import *
from scipy.linalg import eig
from scipy.optimize import minimize

from readsnap import readsnap
from tbplot import *


def main():
    # Set Your Working Directory
    os.chdir('/home/taebum/Workspaces/Intern/MyWork/Datas/m21/')
    # You can set target snapshot number through argv
    # Or inside Program. If argv's length > 2, terminate
    if len(sys.argv) == 1:
        num = int(input("Set Snapshot Number: "))
        until = int(input("How many snapshots to analyze? "))
    elif len(sys.argv) == 2:
        num = int(sys.argv[1])
        until = int(input("How many snapshots to analyze? "))
    elif len(sys.argv) == 3:
        num = int(sys.argv[1])
        until = int(sys.argv[2])
    else:
        exit(1)

    print "Start Snapshot = {:d}".format(num)
    # Only for consecutive snapshots.
    # If you want just one snapshot, input 1
    print "End Snapshot = {:d}\n\n".format(num + until - 1)
    # n is gridsize. Recommanded to be a odd number
    n = 601
    # s_cont: Plot and save contour graph
    # s_fit: Plot and save ellipse fitting compare graph
    # s_angvel: save rotation angle of bar to get angular velocity
    # s_contfit: Plot contour graph and fit graph in one figure
    # s_barfit: Plot bar fitting graph
    s_cont = False
    s_fit = False
    s_append = False
    s_contfit = False
    s_barfit = False
    s_bulge = False
    # If exist bulge at the center, set this parameter true
    e_bulge = False
    if s_append:
        file = open("angle{:d}_{:d}.csv".format(num, num + until - 1), "w")
        writer = csv.writer(file)

    for s_num in range(num, num + until):
        print 'Now in snapshot number {:d}'.format(s_num)

        # Save Snapshot's Data
        total = readsnap(os.getcwd(), s_num, 2, skip_bh=1)
        pos = total['p']
        mass = total['m']
        vel = total['v']

        # Calculate system's center of mass and shift all particles
        tot_mass = sum(mass)
        com = np.array([0, 0, 0], dtype=float)
        com[0], com[1], com[2] = calcom(pos, mass, tot_mass)
        print '\n\nCenter of Mass in x-Domain is {:6.5f}'.format(com[0])
        print 'Mean Position in y-Domain is {:6.5f}'.format(com[1])
        print 'Mean Position in z-Domain is {:6.5f}'.format(com[2])
        print 'Now Shifting These Regions Become 0'
        for i in range(3):
            pos[:, i] -= com[i]

        # Calculate Rho by making grid
        lim = max(abs(min(pos[:, 0])), abs(min(pos[:, 1])), max(pos[:, 0]), max(pos[:, 1])) + 2
        grid = np.mgrid[-lim:lim:complex(0, n)]
        print 'first lim =', lim
        ds = grid[1] - grid[0]
        print '\'ds\' of grid =', ds
        dv = np.power(ds, 3)
        slc = abs(pos[:, 2]) < ds / 2.
        rho = calrho2d(pos[slc, 0], pos[slc, 1], mass[slc], grid, ds)
        rho /= dv
        rho[rho < np.power(np.e, -13)] = np.power(np.e, -13)
        # center density
        rho_bar = rho[n / 2, n / 2]
        # Plot Surface Density Contour
        if s_cont:
            plot_contour(grid, rho, rho_bar, s_num)
        print 'rho_bar(Central Density) =', rho_bar
        # Now Start to Fit ellipse Bar
        # This Fit is Not Exact calculate semi-major and semi-minor length
        # just, to get rotation angle exactly
        # Candidate is arbitrary boundary of ferres bar
        candidate = getcand(grid, rho, rho_bar, 0.25, 0.001)

        e_mat, s21, s22 = setmat(candidate)
        eigen = eig(e_mat, left=False, right=True)
        e_vec = np.array([0, 0, 0], dtype=float)
        for i in range(3):
            e_vec[i] = eigen[1][i][eigen[0] > 0]
        k = 4 * e_vec[2] * e_vec[0] - e_vec[1] * e_vec[1]
        k = np.sqrt(1. / k)
        e_vec *= k
        a3 = -s22 * (s21[0] * e_vec[0] + s21[1] * e_vec[1] + s21[2] * e_vec[2])
        a = np.array([e_vec[0], e_vec[1], e_vec[2], a3])
        fit_chi = alge_dist(candidate[0, :], candidate[1, :], a)
        print 'Ellipse Fitting Finished'
        print 'Ellipse Coefficients: ', a
        print 'Ellipse Fitting chi-square = {:e}'.format(fit_chi)
        theta = 0.5 * np.arctan(a[1] / (a[0] - a[2]))
        major = np.sqrt((-a[3]) / (a[0] * np.power(np.cos(theta), 2) + a[1] * np.cos(theta) * np.sin(theta)
                                   + a[2] * np.power(np.sin(theta), 2)))
        minor = np.sqrt((- a[3]) / (a[0] * np.power(np.sin(theta), 2) - a[1] * np.cos(theta) * np.sin(theta)
                                    + a[2] * np.power(np.cos(theta), 2)))
        if major < minor:
            theta += np.pi / 2
            major, minor = minor, major
        else:
            if theta < 0:
                theta += 2 * np.pi
        print 'theta =', theta, 'rad,', theta * 180. / np.pi, ' degree'
        print 'minor =', minor, 'major =', major

        if s_fit:
            plot_fitcompare(a, candidate, theta, major, minor, fit_chi, s_num)
        if s_append:
            writer.writerow((s_num, theta, theta * 180 / np.pi))
        if s_contfit:
            plot_cont_fit(grid, rho, rho_bar, candidate, a, theta, major, minor, fit_chi, s_num)
        # Rotate Particles to Principal Axes
        pos[:, 0], pos[:, 1] = np.cos(theta) * pos[:, 0] + np.sin(theta) * pos[:, 1], -np.sin(theta) * pos[:, 0] + \
                               np.cos(theta) * pos[:, 1]
        factor = int(lim / 10)
        slc = (abs(pos[:, 0]) < 10.) & (abs(pos[:, 1]) < 10.) & (abs(pos[:, 2]) < 10.)
        lim = 10
        n /= factor
        n += 1
        # lim = max(abs(min(pos[:, 0])), abs(min(pos[:, 1])), max(pos[:, 0]), max(pos[:, 1]))
        grid = np.mgrid[-lim:lim:complex(0, n)]
        ds = grid[1] - grid[0]
        print '\'ds\' of grid =', ds
        dv = np.power(ds, 3)
        print 'second boxsize =', n
        print 'second lim =', lim
        rho, index = calrho3d(pos[slc], mass[slc], grid, ds, 1)
        rho /= dv
        rho_bar = rho[n / 2][n / 2][n / 2]
        print 'rho bar in 3-dim =', rho_bar

        rho[np.where(rho < np.power(np.e, -13))] = np.power(np.e, -13)
        if e_bulge:
            print 'Bulge Existing'
            bulge = np.array([1], dtype=float)
            bulge_result = minimize(compbulge, bulge, args=(pos[slc], index, rho, rho_bar), method='SLSQP')
            bulge = bulge_result['x']
            print bulge_result
            if s_bulge:
                plot_bulge_fit(grid, rho, rho_bar, bulge[0], s_num)

        vector = np.array([1, major, minor])
        result = minimize(compbar, vector, args=(pos[slc], index, rho, rho_bar), method='SLSQP')
        print result
        vector = result['x']

        newvector = np.array([rho_bar * 0.7, 0.5, 0.7, 4, 1])
        print newvector
        newresult = fitbulge(newvector, pos[slc], index, rho, rho_bar)
        print newresult

        mimi = minimize(fitbulge, newvector, args=(pos[slc], index, rho, rho_bar), method='SLSQP')
        print mimi
        newvector = mimi['x']
        x_sli = rho[n / 2, n / 2, :]
        y_sli = rho[n / 2, :, n / 2]
        z_sli = rho[:, n / 2, n / 2]
        amaj = np.linspace(-newvector[3], newvector[3], 100)
        amin = np.linspace(-newvector[4], newvector[4], 100)
        rho_bulge = newvector[0]
        rho_ferres = rho_bar - newvector[0]
        b = np.linspace(-8, 8, 100)
        fig, ax = plt.subplots(1, 3, sharey='row', figsize=(16, 9))
        ax[0].plot(grid, x_sli)
        ax[0].plot(b, rho_bulge * np.power(1 + np.power(b / newvector[1], 2), -1.5))
        ax[0].plot(amaj, rho_ferres * np.power(1 - np.power(amaj / newvector[3], 2), newvector[2]))
        ax[0].set_xlim(-8, 8)
        ax[0].set_xticks(np.mgrid[-8:8:5j])
        ax[0].set_ylim(0, 1.2)
        ax[0].tick_params(axis='y', labelsize=15)
        ax[0].tick_params(axis='x', labelsize=15)
        text = 'Snapshot = {:d}\n'.format(s_num)
        text += r'$\rho_{bul}$ = ' + '{:5.4f}, '.format(newvector[0])
        text += r'$\rho_{bar}$ = ' + '{:5.4f}\n'.format(rho_bar - newvector[0])
        text += r'$R_{bul}$ = ' + '{:5.4f} '.format(newvector[1])
        text += 'n = {:5.4f}\n'.format(newvector[2])
        text += 'major = {:5.4f}\n'.format(newvector[3])
        text += 'major = {:5.4f}'.format(newvector[4])
        ax[0].text(-7.8, 1.2 * 0.8, text, fontsize=18)
        ax[1].plot(grid, y_sli)
        ax[1].plot(b, rho_bulge * np.power(1 + np.power(b / newvector[1], 2), -1.5))
        ax[1].plot(amin, rho_ferres * np.power(1 - np.power(amin / newvector[4], 2), newvector[2]))
        ax[1].tick_params(axis='x', labelsize=15)
        ax[1].set_xlim(-8, 8)
        ax[1].set_xticks(np.mgrid[-8:8:5j])
        ax[2].plot(grid, z_sli, label='Real(Summation)')
        ax[2].plot(b, rho_bulge * np.power(1 + np.power(b / newvector[1], 2), -1.5), label='Bulge')
        ax[2].plot(amin, rho_ferres * np.power(1 - np.power(amin / newvector[4], 2), newvector[2]), label='Bar')
        ax[2].tick_params(axis='x', labelsize=15)
        ax[2].set_xlim(-8, 8)
        ax[2].set_xticks(np.mgrid[-8:8:5j])
        ax[2].legend(fontsize=15)
        fig.tight_layout()
        plt.savefig('barbar_{:d}.png'.format(s_num))
        plt.close()

        if s_barfit:
            x_sli = rho[n / 2, n / 2, :]
            y_sli = rho[n / 2, :, n / 2]
            z_sli = rho[:, n / 2, n / 2]
            a = np.linspace(-vector[1] + 1e-6, vector[1] - 1e-6, 100)
            b = np.linspace(-vector[2] + 1e-6, vector[2] - 1e-6, 100)

            plot_bar_fit(grid, x_sli, y_sli, z_sli, rho_bar, a, b, vector, s_num)


def getdist(x, y, z):
    return np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))


def alge_dist(xx, yy, u):
    result = 0
    amount = len(xx)
    for i in range(amount):
        result += np.power(u[0] * np.power(xx[i], 2) + u[1] * xx[i] * yy[i] + u[2] * np.power(yy[i], 2) + u[3], 2)
    return result / amount


def getcand(grid, rho, rho_bar, percentage, eps):
    plt.figure()
    cc = plt.contour(grid, grid, rho, levels=np.mgrid[rho_bar * (percentage - eps):rho_bar * (percentage + eps):256j])
    plt.gcf().clear()
    plt.close()
    length = 0
    where = 0
    # Find The Most Particles Level
    for i in range(len(cc.collections[0].get_paths())):
        if len(cc.collections[0].get_paths()[i]) > length:
            length = len(cc.collections[0].get_paths()[i])
            where = i
    pp = cc.collections[0].get_paths()[where]
    vv = pp.vertices
    return np.array((vv[:, 0], vv[:, 1]), dtype=float)


if __name__ == "__main__":
    main()
