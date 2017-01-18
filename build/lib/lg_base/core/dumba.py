
def _isFunc(maybefunc):
	t=type(maybefunc)
	#hasattr(args[0],'__call__')
	return t==type(_isFunc) or t==type or t==type(doNothingDecorator)

def returnParm(parm):
	return parm

class aNumbaType(object):
	def __init__(self):
		pass

	def __call__(self,*args,**kwargs):
		return returnParm #in the case its being used as a decorator for an object's function signature

	def __getitem__(self,idx):
		return self

#a decorator that is always called
def alwaysCalledDoNothingDecorator(*l1args,**l1kwargs):
	return returnParm

#a decorator that is never called
def neverCalledDoNothingDecorator(*l1args,**l1kwargs):
	if len(l1kwargs)!=0 or len(l1args)!=1:
		#definately a called decorator
		raise RuntimeError('Called a nevercalled')
	return l1args[0]

#a decorator that is either not called, or called in a way that is obvious
def neverSimplyCalledDoNothingDecorator(*l1args,**l1kwargs):
	if len(l1kwargs)!=0 or len(l1args)!=1 or isinstance(l1args[0],aNumbaType) or l1args[0]==returnParm:
		#definately a called decorator
		return returnParm
	return l1args[0]#simple 1 parameter, so not called, this is the function

def doNothingDecorator(*l1args,**l1kwargs):
	def l2decorator(*l2args,**l2kwargs):
		if len(l2kwargs)!=0 or len(l2args)!=1:
			#definately a called here. decorator couldn't have been called
			return l1args[0](*l2args,**l2kwargs)
		if _isFunc(l2args[0]):
			return l2args[0]
		#probably not called. just return first param
		print 'Warning: NOT SURE WHAT TO DO in doNothingDecorator'
		raise RuntimeError("Decorator couldn't figure out if it was called or not")
		return l1args[0]

	if len(l1kwargs)!=0 or len(l1args)!=1:
		#definately a called decorator
		return returnParm
	#only one parameter, and not named. could be a function
	datFunc=l1args[0]
	if isinstance(datFunc,basestring):#called with a string
		return returnParm
	if datFunc==alwaysCalledDoNothingDecorator or isinstance(datFunc,aNumbaType) or datFunc==returnParm:
		return returnParm
	if _isFunc(datFunc):
		return datFunc
	return l2decorator

void=aNumbaType()
uint16=aNumbaType()
int_=aNumbaType()
float_=aNumbaType()
float64=aNumbaType()
double=aNumbaType()
object_=aNumbaType()
autojit=neverSimplyCalledDoNothingDecorator
jit=doNothingDecorator
generated_jit=doNothingDecorator
jitclass=doNothingDecorator
vectorize=doNothingDecorator
