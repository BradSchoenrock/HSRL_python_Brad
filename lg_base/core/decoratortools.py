import os
from datetime import datetime

def package_version(name):
	import pkg_resources
	try:
		v=pkg_resources.require(name)[0]
	except pkg_resources.DistributionNotFound:
		return None
	return v.version

def package_min_version(name,minvers):
	minvarr=minvers.split('.')
	v=package_version(name)
	if v is None:
		return False
	versarr=v.split('.')
	for i,x in enumerate(minvarr):
		if x>versarr[i]:
			return False
		elif x<versarr[i]:
			break
		else: #same. on to next
			pass
	return True

def package_max_version(name,maxvers):
	maxvers=maxvers.split('.')
	v=package_version(name)
	if v is None:
		return False
	versarr=v.split('.')
	for i,x in enumerate(maxvers):
		if x<versarr[i]:
			return False
		elif x>versarr[i]:
			break
		else: #same. on to next
			pass
	return True

class Delayed(object):
	def __init__(self,decorator,minversion=None,maxversion=None,runAtImport=False,timeDecoration=False):
		self.dec=decorator
		self.realfunc=None
		self.canRun=True
		self.runAtImport=runAtImport
		self.timeDecoration=timeDecoration
		if minversion!=None:
			for pkgname,minvers in minversion.items():
				if not package_min_version(pkgname,minvers):
					self.canRun=False
					print 'version of',pkgname,'is',package_version(pkgname),'which is less than',minvers

		if maxversion!=None:
			for pkgname,maxvers in maxversion.items():
				if not package_max_version(pkgname,maxvers):
					self.canRun=False
					print 'version of',pkgname,'is',package_version(pkgname),'which is greater than',maxvers


	def calldecorator(self):
		if self.timeDecoration:
			starttime=datetime.utcnow()
			print 'running delayed decorator',self.dec,'on',self.func
		self.realfunc=self.dec(self.func)
		if self.timeDecoration:
			print 'took',(datetime.utcnow()-starttime).total_seconds(),'seconds'

	def docall(self,*args,**kwargs):
		if self.realfunc==None:
			self.calldecorator()
		return self.realfunc(*args,**kwargs)

	def __call__(self,func):
		if not self.canRun:
			print 'Dropping decorator',self,'for',func
			return func
		self.func=func
		if self.runAtImport or os.getenv('NODELAYDECORATOR',None)!=None:
			self.calldecorator()
			return self.realfunc
		return self.docall
