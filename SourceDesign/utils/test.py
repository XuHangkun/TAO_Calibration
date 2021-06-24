
class Parent:
    def hello(self):
        print('我是父亲')


class Child(Parent):
    # child类继承与 parent类
    def hello(self):
        # parent类中有hello方法，但是这里也定义了一个hello方法。
        # 覆盖,但是父类方法不受影响
        print('我是孩子')


def main():
    Child().hello()  # 调用的是 子类的方法
    Parent().hello()  # 调用的是 父类的方法


if __name__ == '__main__':
    main()