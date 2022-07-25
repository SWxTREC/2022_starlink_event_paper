import matplotlib.pyplot as plt
import numpy as np

# Starlink launches retrieved from:
# https://en.wikipedia.org/wiki/List_of_Starlink_launches

# Starlink obtained from space-track.org
# https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/51460/orderby/TLE_LINE1 ASC/EPOCH/2022-02-03--2022-02-11/format/tle
# http: // www.celestrak.com/Norad/elements/table.php?tleFile = starlink & title = Starlink % 20Satellites & orbits = 0 & pointsPerRev = 90 & frame = 1

# F107 and Ap data retrieved from:
# https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt
data = """2019-05-24T02:30,68.1,71.4,4.0,2.6e-11
2019-11-11T14:56,69.3,67.9,3.0,2.8e-11
2020-02-17T15:05,69.4,69.6,13.0,2.9e-11
2020-03-18T12:16,71.5,69.6,10.0,3.1e-11
2020-04-22T19:30,69.5,70.3,3.0,2.8e-11
2020-06-04T01:25,72.1,71.0,4.0,2.6e-11
2020-06-13T09:21,71.6,71.2,2.0,2.4e-11
2020-08-07T05:12,76.0,72.1,3.0,2.5e-11
2020-08-18T14:31,72.2,72.3,5.0,2.6e-11
2020-09-03T12:46,70.8,72.6,7.0,2.8e-11
2020-10-06T11:29,71.6,74.9,6.0,3.1e-11
2020-10-18T12:25,74.2,77.8,6.0,3.3e-11
2020-10-24T15:31,73.3,80.2,14.0,3.5e-11
2020-11-25T02:13,101.0,83.9,6.0,3.6e-11
2021-01-20T13:02,75.2,74.8,1.0,2.6e-11
2021-01-24T15:00,74.8,74.5,15.0,3e-11
2021-02-04T06:19,71.2,73.9,7.0,2.8e-11
2021-02-16T03:59,69.1,73.2,12.0,3e-11
2021-03-04T08:24,74.1,73.1,9.0,3.1e-11
2021-03-11T08:13,74.6,73.4,3.0,2.9e-11
2021-03-14T10:01,74.2,73.6,23.0,3.4e-11
2021-03-24T08:28,75.3,74.3,9.0,3.2e-11
2021-04-07T16:34,70.5,74.8,4.0,3e-11
2021-04-29T03:44,78.0,76.2,4.0,3.1e-11
2021-05-04T19:01,70.7,76.4,1.0,2.8e-11
2021-05-09T06:42,75.0,76.8,3.0,2.9e-11
2021-05-15T22:56,74.9,77.6,4.0,2.9e-11
2021-05-26T18:59,85.5,80.4,14.0,3.3e-11
2021-06-30T19:31,97.2,81.8,5.0,2.9e-11
2021-09-14T03:55,79.0,84.2,5.0,3.2e-11
2021-11-13T11:19,81.9,89.3,1.0,3.3e-11
2021-12-02T23:12,82.9,93.8,8.0,3.6e-11
2021-12-18T12:41,111.6,94.6,14.0,4e-11
2022-01-06T21:49,103.8,102.8,0.0,3.5e-11
2022-01-19T02:02,101.9,106.5,23.0,4.3e-11
2022-02-03T18:13,126.0,106.4,32.0,4.8e-11"""

dates = []
F107s = []
F107as = []
Aps = []
densities = []

for line in data.split("\n"):
    date, F107, F107a, Ap, density = line.split(",")
    dates.append(np.datetime64(date))
    F107s.append(float(F107))
    F107as.append(float(F107a))
    Aps.append(float(Ap))
    densities.append(float(density))

plt.rcParams.update({"text.usetex": True, "font.family": "serif"})
fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, sharex=True,
                                    constrained_layout=True)


densities = [d/1e-11 for d in densities]
ax1.plot(dates, densities, 'o-.', c='k')
ax1.set_yticks([2, 3, 4, 5])
ax1.axhline(np.mean(densities[:-1]), c='k', alpha=0.3)
ax1.set_ylabel("Density ($10^{-11}$ kg/m$^3$)")
ax2.plot(dates, F107s, 'o-.', c='k')
ax2.axhline(np.mean(F107s[:-1]), c='k', alpha=0.3)
ax2.set_yticks([75, 100, 125])
ax2.set_ylabel("F10.7 (sfu)")
ax3.plot(dates, Aps, 'o-.', c='k')
ax3.axhline(np.mean(Aps[:-1]), c='k', alpha=0.3)
ax3.set_yticks([0, 10, 20, 30, 40])
ax3.set_ylabel("Ap")
fig.suptitle("Space weather conditions for Starlink launches")
fig.align_labels()
# plt.show()
fig.savefig("starlink-conditions.pdf")
