import time

class GlobalTimer:
    _start = None
    
    @classmethod
    def start(cls):
        cls._start = time.time()
    
    @classmethod
    def elapsed(cls):
        if cls._start is None:
            cls.start()
        return time.time() - cls._start
    
    @classmethod
    def log(cls, message):
        elapsed = cls.elapsed()
        print(f"[{elapsed:7.2f}s] {message}")
