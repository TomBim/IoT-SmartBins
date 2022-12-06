import sys
sys.platform
import rpyc

c = rpyc.classic.connect("c2")
c.modules.sys