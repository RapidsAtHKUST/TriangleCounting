import datetime
import os

if __name__ == '__main__':
    date_str = datetime.datetime.now().strftime("%Y-%m-%d-%H:%M")
    print(date_str)
    os.system('zip -r tc-rapids-{}.zip triangle-counting -x cmake-build-debug -x */CMake*'.format(date_str))
