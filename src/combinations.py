# taken from Python cookbook, 2nd Edition, Page 725 - http://books.google.com/books?id=Q0s6Vgb98CQC&pg=PT759&lpg=PT759&dq=Python+cookbook+_combinators&source=bl&ots=hc3-5aPmtG&sig=X6_DsalZqYwn7zA_0vMhlj9lqmU&hl=en&ei=OtZqS6SWDIvQtAOii4CgAw&sa=X&oi=book_result&ct=result&resnum=3&ved=0CBIQ6AEwAg#v=onepage&q=&f=false

def _combinators(_handle, items, n):
    if n==0:
        yield []
        return
    for i, item in enumerate(items):
        this_one = [ item ]
        for cc in _combinators(_handle, _handle(items, i), n-1):
            yield this_one + cc

def selections(items, n):
    # take n (not necessarily distinct) items, order matters
    def keepAllItems(items, i):
        return items
    return _combinators(keepAllItems, items, n)

from itertools import product, izip, tee

def iproduct(items, n):
    try:
        return product(items, repeat = n)
    # if python < 2.6
    except AttributeError:
        return selections(items, n) # from src.combinations.py (Python cookbook)

# from: http://www.python.org/doc//current/library/itertools.html#recipes
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)
