import matplotlib.pyplot as plt
import graphs

x = [1.60944, 2.30259, 2.99573, 3.68888, 4.38203, 5.07517]
y1 = [8.36058, 7.07201, 5.67014, 4.25532, 2.849, 1.43123]
y2 = [8.48089, 6.42638, 4.73388, 3.23505, 1.80579, 0.388583]
fig, ax = graphs.basePlot()
plt.title("log(error) over log(N)")
plt.xlabel("log(N)")
plt.ylabel("log(error)")
ax.plot(x, y1, label = "Natural: k ~ -2")
ax.plot(x, y2, label = "Unnatural: 1й уч-к (до 2.3): k ~ -3; 2й: k ~ -2")
plt.legend()
plt.savefig('../CubicSpline_log(error)(log(N)).png')
plt.show()