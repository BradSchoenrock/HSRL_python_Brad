import traceback
import time
import os,datetime

def findGit(path):
    if not os.path.exists(path):
        return None
    path=os.path.realpath(path)
    oldpath=None
    while oldpath!=path:
        if os.path.exists(os.path.join(path,'.git')):
            return path
        oldpath=path
        path=os.path.dirname(path)
    return None

def getCurrentHash(path):
    path=findGit(path)
    if path is None:
        return None,None
    try:
        import git
    except ImportError:
        return None,None
    try:
        rep=git.Repo(path)
    except git.InvalidGitRepositoryError:
        raise ValueError('Not a GIT tracked location')
    except:
        print 'Unknown API for git'
        traceback.print_exc()
        raise RuntimeError('Unknown Git command. FIXME')

    try:
        commit=rep.commits('HEAD')[0]
        commit_date=datetime.datetime(*(commit.committed_date[:6]))
        commitid=commit.id
    except AttributeError: #is 0.3 or newer
        try:
            commitlist=rep.iter_commits('HEAD')
            print commitlist
            commit=None
            for c in commitlist:
                if c!=None:
                    commit=c
                    break
            if commit==None:
                return None,None
            commit_date=datetime.datetime(*time.gmtime(commit.committed_date)[:6])
            commitid=commit.name_rev.split(' ')[0]
        except git.BadObject:
            return None,None
        except:
            print 'Unknown API for git'
            traceback.print_exc()
            raise RuntimeError('Unknown Git command. FIXME')
    
    return commitid,commit_date

def getCodeVersion():
    codevers=None
    codedate=None
    try:
        codevers,codedate=getCurrentHash(os.path.dirname(__file__))
    except ValueError:
        pass
    if codevers==None:#no git version
        import pkg_resources
        try:
            v=pkg_resources.require('hsrl')[0]
            codevers=v.version
        except pkg_resources.DistributionNotFound as e:
            pass
    if codevers==None:
        codevers="unknown"
    return codevers,codedate

def main():
    import sys
    if len(sys.argv)==1:
        print getCodeVersion()
    else:
        for path in sys.argv[1:]:
            print path
            print getCurrentHash(path)

if __name__ == '__main__':
    main()
