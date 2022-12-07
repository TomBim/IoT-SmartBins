import rpyc
from rpyc.utils.server import ThreadedServer


@rpyc.service
class MyService(rpyc.Service):
    print("EITA PEGAAAAA")
    @rpyc.exposed
    def add(self,a,b):
        print(a+b)
        return a+b
    
    @rpyc.exposed
    def sub(self,a,b):
        return a-b

    @rpyc.exposed
    def mult(self,a,b):
        return a*b
    
    @rpyc.exposed
    def div(self,a,b):
        if b==0:
            print("safadooo... n pode dividir por zero, tá de sacanagem???")
        else:
            return a/b
    
    def foo(self):
        print("foo")


if __name__ == "__main__":
    server = ThreadedServer(MyService, port=12345)
    server.start()