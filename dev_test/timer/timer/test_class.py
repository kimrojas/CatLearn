import time


class Timer:
    def __init__(self):
        self.start_time = None
        self.stop_time = None

    def start(self, name):
        self.start_time = time.time()
        print(f"Timer '{name}' started")

    def stop(self, name):
        self.stop_time = time.time()
        print(f"Timer '{name}' stopped")

    def timer_decorator(self, name):
        def decorator_func(func):
            def wrapper(*args, **kwargs):
                self.start(name)
                result = func(*args, **kwargs)
                self.stop(name)
                return result

            return wrapper

        return decorator_func

    def report(self):
        if self.start_time and self.stop_time:
            elapsed_time = self.stop_time - self.start_time
            print(f"Elapsed time: {elapsed_time:.6f} seconds")
        else:
            print("Timer has not been started or stopped yet")


# Base class for MyClass
class MyBaseClass:
    def base_method(self):
        # Base class method code
        pass


# MyClass with Timer functionality
class MyClass(MyBaseClass):
    def __init__(self):
        super().__init__()
        self.timer = Timer()  # Create a Timer instance

    @self.timer.timer_decorator(name="my_method")
    def my_method(self):
        # Your method code here
        time.sleep(2)


# Usage example
