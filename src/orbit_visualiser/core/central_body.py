

class CentralBody():

    def __init__(self):
        self._mu = 398600 # gravitational parameter in km³/s² = Gm
        self._r = 6378 # radius in km

    @property
    def mu(self) -> float:
        return self._mu

    @mu.setter
    def mu(self, mu: float) -> None:
        self._mu = mu