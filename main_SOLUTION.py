import sys
sys.platform
import rpyc

conn = rpyc.connect("0.0.0.0", 12345)
x = conn.root.add(4,7)
y
print(x)