import matplotlib.pyplot as plt
import graphs

x = [2, 1, 0.5, 0.25, 0.125, 0.0625]
y1 = [349.089, 39.1576, 6.57764, 1.34839, 0.305257, 0.07262]
y2 = [101.331, 31.0841, 12.3521, 5.52109, 2.6116, 1.27025]
y3 = [222.827, 23.6963, 6.04527, 2.32237, 1.00661, 0.467012]
fig, ax = graphs.basePlot()
plt.title("error over length")
plt.xlabel("length")
plt.ylabel("error")
ax.plot(x, y1, label = "N = 3")
ax.plot(x, y2, label = "N = 4")
ax.plot(x, y3, label = "N = 5")
plt.legend()
plt.savefig('../HermiteInterpolator_error(length).png')
plt.show()