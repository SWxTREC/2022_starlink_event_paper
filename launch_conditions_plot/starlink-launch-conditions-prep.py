import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
from pymsis import msis
from skyfield.api import EarthSatellite, load, wgs84

# Starlink launches
# https://en.wikipedia.org/wiki/List_of_Starlink_launches
launch_dates = ["2019-05-24T02:30",
                "2019-11-11T14:56",
                "2020-02-17T15:05",
                "2020-03-18T12:16",
                "2020-04-22T19:30",
                "2020-06-04T01:25",
                "2020-06-13T09:21",
                "2020-08-07T05:12",
                "2020-08-18T14:31",
                "2020-09-03T12:46",
                "2020-10-06T11:29",
                "2020-10-18T12:25",
                "2020-10-24T15:31",
                "2020-11-25T02:13",
                "2021-01-20T13:02",
                "2021-01-24T15:00",
                "2021-02-04T06:19",
                "2021-02-16T03:59",
                "2021-03-04T08:24",
                "2021-03-11T08:13",
                "2021-03-14T10:01",
                "2021-03-24T08:28",
                "2021-04-07T16:34",
                "2021-04-29T03:44",
                "2021-05-04T19:01",
                "2021-05-09T06:42",
                "2021-05-15T22:56",
                "2021-05-26T18:59",
                "2021-06-30T19:31",
                "2021-09-14T03:55",
                "2021-11-13T11:19",
                "2021-12-02T23:12",
                "2021-12-18T12:41",
                "2022-01-06T21:49",
                "2022-01-19T02:02",
                "2022-02-03T18:13",
                ]

text = """
STARLINK-3167           
1 51460U 22010E   22042.23745443  .00315859  49990-4  77952-3 0  9992
2 51460  53.2159 135.9152 0033271 206.0443 153.8913 15.94471269  2420
"""

# Starlink obtained from space-track.org
# https://www.space-track.org/basicspacedata/query/class/gp_history/NORAD_CAT_ID/51460/orderby/TLE_LINE1 ASC/EPOCH/2022-02-03--2022-02-11/format/tle
# http: // www.celestrak.com/Norad/elements/table.php?tleFile = starlink & title = Starlink % 20Satellites & orbits = 0 & pointsPerRev = 90 & frame = 1
text = """
STARLINK
1 51460U 22010E   22039.75001157 -.00949494  29636-3 -17121-2 0  9997
2 51460  53.2172 148.6578 0053662 179.7290 281.5305 15.98078309  2033
"""
lines = text.strip().splitlines()

sat = EarthSatellite(lines[1], lines[2], lines[0])
ts = load.timescale()
daily_times = ts.utc(2022, 2, 5, 0, range(1440))

t = ts.now()
# geocentric = sat.at(t)
geocentric = sat.at(daily_times)

pos = wgs84.geographic_position_of(geocentric)

lons = pos.longitude.degrees
lats = pos.latitude.degrees
alts = pos.elevation.km

# F107 and Ap data
# https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt
dates = daily_times.utc_datetime()
ndates = len(dates)
f107 = np.full(ndates, 80)
f107a = np.full(ndates, 80)
aps = np.full((ndates, 7), 4)

storm_output = np.ones((ndates, 11))
f107_only = np.ones_like(storm_output)
ap_only = np.ones_like(storm_output)
quiet_output = np.ones_like(storm_output)
previous_launch = np.ones_like(storm_output)

for i in range(ndates):
    lon = lons[i]
    lat = lats[i]
    alt = alts[i]
    # Previous launch (2022-01-19)
    previous_launch[i, :] = msis.run(dates[i], lon, lat, alt, 105, 105, [[23]*7])
    # Quiet conditions
    quiet_output[i, :] = msis.run(dates[i], lon, lat, alt, 80, 80, [[4]*7])
    # Increase f107
    f107_only[i, :] = msis.run(dates[i], lon, lat, alt, 125, 125, [[4]*7])
    # Increase Ap
    ap_only[i, :] = msis.run(dates[i], lon, lat, alt, 80, 80, [[36]*7])
    # Storm conditions (f107 and Ap)
    storm_output[i, :] = msis.run(dates[i], lon, lat, alt, 125, 125, [[36]*7])

fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, constrained_layout=True)
plt.rcParams.update({"text.usetex": True, "font.family": "serif"})
import matplotlib.units as munits
import matplotlib.dates as mdates
import datetime
converter = mdates.ConciseDateConverter()
munits.registry[np.datetime64] = converter
munits.registry[datetime.date] = converter
munits.registry[datetime.datetime] = converter

densities = [d/1e-11 for d in previous_launch[:, 0]]
ax1.plot(dates, densities, c='k')
ax1.set_xlim(dates[0], dates[-1])
ax1.set_ylabel("Density ($10^{-11}$ kg/m$^3$)")
ax2.plot(dates, alts, c='k')
ax2.set_ylabel("Altitude (km)")
fig.suptitle("Initial orbit of Starlink satellite")
fig.align_labels()
# ax1.set_yticks([2, 3, 4, 5])
fig.savefig("starlink-altitude.pdf")

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


previous_launch_data = moving_average(previous_launch[:, 0], 90)
quiet_data = moving_average(quiet_output[:, 0], 90)
storm_data = moving_average(storm_output[:, 0], 90)
f107_data = moving_average(f107_only[:, 0], 90)
ap_data = moving_average(ap_only[:, 0], 90)

# fig, ax = plt.subplots()
# ax.plot(dates[45:-44], quiet_data, label="F107=80, Ap=4")
# ax.plot(dates[45:-44], f107_data, label="F107=125, Ap=4")
# ax.plot(dates[45:-44], ap_data, label="F107=80, Ap=36")
# ax.plot(dates[45:-44], storm_data, label="F107=125, Ap=36")
# ax.plot(dates[45:-44], previous_launch_data, c="k", label="Jan 19th 2022")
# ax.plot(dates[45:-44], quiet_data/quiet_data, label="F107=80, Ap=4")
# ax.plot(dates[45:-44], f107_data/quiet_data, label="F107=125, Ap=4")
# ax.plot(dates[45:-44], ap_data/quiet_data, label="F107=80, Ap=36")
# ax.plot(dates[45:-44], storm_data/quiet_data, label="F107=125, Ap=36")
# ax.plot(dates[45:-44], previous_launch_data/quiet_data,
#         c="k", label="Jan 19th 2022")
# ax.legend()
# plt.show()

# ax.plot(dates, quiet_output[:, 0], label="F107=80, Ap=4")
# ax.plot(dates, f107_only[:, 0] / quiet_output[:, 0], label="F107=125, Ap=4")
# ax.plot(dates, ap_only[:, 0] / quiet_output[:, 0], label="F107=80, Ap=36")
# ax.plot(dates, storm_output[:, 0] /
#         quiet_output[:, 0], label="F107=125, Ap=36")
# ax.plot(dates, previous_launch[:, 0] /
#         quiet_output[:, 0], c="k", label="Jan 19th 2022")
# ax.legend()
# plt.show()

quiet_mean = quiet_output[:, 0].mean()
f107_mean = f107_only[:, 0].mean()
ap_mean = ap_only[:, 0].mean()
storm_mean = storm_output[:, 0].mean()
previous_launch_mean = previous_launch[:, 0].mean()

print("Differences:")
print(f"F107: {f107_mean / quiet_mean}")
print(f"Ap: {ap_mean / quiet_mean}")
print(f"Storm: {storm_mean / quiet_mean}")
print(f"Previous launch: {previous_launch_mean / quiet_mean}")



def xyz(lon, lat, r):
    """Convert spherical to Cartesian coordinates."""
    r += 6378
    x = r * np.cos(np.deg2rad(lon)) * np.cos(np.deg2rad(lat))
    y = r * np.sin(np.deg2rad(lon)) * np.cos(np.deg2rad(lat))
    z = r * np.sin(np.deg2rad(lat))
    return x, y, z

x, y, z = xyz(lons, lats, alts)

ax = plt.axes(projection='3d')

ax.plot(x, y, z)

plt.show()

# Actual data

df = pd.read_csv("Kp_ap_Ap_SN_F107_since_1932.txt", 
                 delim_whitespace=True,
                 parse_dates={"date": ["year", "month", "day"]},
                 na_values=["-1", "-1.0", "-1.000"],
                 skiprows=40,
                 header=None,
                 index_col="date",
                 usecols=[0, 1, 2, 23, 26],
                 names=["year", "month", "day", "Ap", "F10.7"])
# There are some anomalous 0s in F10.7 as well, so get rid of them too
df.loc[df["F10.7"] < 1, "F10.7"] = np.nan

df["F10.7a"] = df["F10.7"].rolling(81, center=True, min_periods=1).mean()


def get_f107_data(date):
    """Get the f107 and ap data closest to the given date"""
    closest_index = np.argmin(np.abs(df.index - date))
    F107 = df["F10.7"][closest_index]
    F107a = df["F10.7a"][closest_index]
    Ap = df["Ap"][closest_index]
    return (F107, F107a, Ap)

nevents = len(launch_dates)
event_densities = np.zeros(nevents)

print(f"date,F107,F107a,Ap")
for j, date in enumerate(launch_dates):
    f107, f107a, Ap = get_f107_data(np.datetime64(date))
    # Temporary fill values
    if np.isnan(f107):
        f107 = 75

    event_times = ts.utc(int(date[:4]), int(date[5:7]), int(
        date[8:10]), 0, range(1440)).utc_datetime()
    ndates = len(event_times)
    output = np.ones((ndates, 11))
    for i in range(ndates):
        lon = lons[i]
        lat = lats[i]
        alt = alts[i]
        # Previous launch (2022-01-19)
        output[i, :] = msis.run(
            event_times[i], lon, lat, alt, f107, f107a, [[Ap]*7])

    event_mean = output[:, 0].mean()
    event_densities[j] = event_mean
    print(f"{date},{f107:.1f},{f107a:.1f},{Ap},{event_mean:.2g}")

fig, ax = plt.subplots()

ax.plot(launch_dates, event_densities)

plt.show()
