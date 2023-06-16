from catlearn.utilities.timings import Timer
from time import time, sleep

t = Timer(debug=True, logfile="test_timer.log")
print(t.get_time())

t.start("test start")
sleep(2)
t.stop("test start")

t.report()


@t.timer_decorator(name="test_func")
def test_func():
    sleep(1)


test_func()
t.report()
