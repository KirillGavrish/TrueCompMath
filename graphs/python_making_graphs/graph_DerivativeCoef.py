import matplotlib.pyplot as plt
import graphs

x = [0, -2.30259, -4.60517, -6.90776, -9.21034, -11.5129, -13.8155, -16.1181, -18.4207, -20.7233, -23.0259, -25.3284, -27.631, -29.9336, -32.2362, -34.5388]
y1 = [1.66911, -1.24345, -3.59933, -5.90717, -8.21037, -10.8548, -7.74636, -2.3212, 1.81865, 6.10213, 1, 15.3064, 20.6047, 24.5167, 29.1219, 33.727]
y2 = [2.13806, -3.58179, -8.28643, -12.9018, -16.3139, -10.2898, -7.74636, -1.67679, 2.71106, 6.10213, 10.7013, 1, 22.1088, 25.2099, 29.1219, 35.1133]
y3 = [2.58311, -5.92966, -12.9819, -19.0392, -17.0908, -10.2898, -7.74636, -0.612038, 2.36127, 8.17623, 12.0875, 17.7913, 22.5506, 24.5167, 29.1219, 34.4202]
fig, ax = graphs.basePlot()
plt.title("log(delta) over log(h)")
plt.xlabel("log(h)")
plt.ylabel("log(delta)")
ax.plot(x, y1, label = "N = 3, 1й уч-к: k ~ -2, 2й: k ~ 1.14")
ax.plot(x, y2, label = "N = 4, 1й уч-к: k ~ -2, 2й: k ~ 2.26")
ax.plot(x, y3, label = "N = 5, 1й уч-к: k ~ -2, 2й: k ~ 3.38")
plt.legend()
plt.savefig('../DerivativeCoef_log(delta)(log(h)).png')
plt.show()
