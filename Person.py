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
    _map_=None

    def __init__(self, id, origin, destination, map_, *args):
        self._id = id
        self._origin = origin
        self._destination = destination
        self._map_ = map_
        self._path = NotImplemented

    def _path_planner(self):
        pass

    def get_id(self):
        return self._id

    def get_path(self):
        return self._path

    def get_origin(self):
        return self._origin

    def get_destination(self):
        return self._destination
    
    def set_origin(self, origin):
        self._origin = origin

    def set_destination(self, destination):
        self._destination = destination
