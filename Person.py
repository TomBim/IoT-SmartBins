from ssl import OP_NO_RENEGOTIATION


class Person:
    _id=None
    _pos=None
    _origin=None
    _destination=None
    _path=None #list of intersection points
    _fov=None
    _max_distance_carrying_trash=None #maximal distance that this person carries a trash
    _time_of_consumption=None #time needed for consuming
    _speed=None
    _has_trash=False #if consuming something, it's false
    _time_alive=None
    _distance_carrying_trash=None

    def __init__(self, id, origin, destination, education):
        self._id = id
        self._origin = origin
        self._destination = destination
        self._education = education
        self._path = ??? #create path

    def getID(self):
        return self._id

    def getEducation(self):
        return self._education

    def getPath(self):
        return self._path

    def getOrigin(self):
        return self._origin

    def getDestination(self):
        return self._destination
    
    def setOrigin(self, origin):
        self._origin = origin
        self._path = ???

    def setDestination(self, destination):
        self._destination = destination
        self._path = ???
