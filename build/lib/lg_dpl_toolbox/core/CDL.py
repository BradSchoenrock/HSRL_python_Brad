
from ply import lex,yacc
import sys,os
from collections import OrderedDict
import numpy

class CDLLexer:
    tokens = (
        'LETTERS',
        'DIGITS',
        'SIGN',
        'PERIOD',
        'SEMICOLON',
        'EQUALS',
        'OPENPAREN',
        'CLOSEPAREN',
        'OPENBRACE',
        'CLOSEBRACE',
        'OPENBRACKET',
        'CLOSEBRACKET',
        'COLON',
        'QUOTE',
        'COMMA',
        'WHITESPACE',
        'UNDERSCORE',
        'SLASH',
        'BACKSLASH',
        'LINEEND',
        'OTHER',
        )

    def t_LINEEND(self,t):
        r'(\n|\r|\n\r|\r\n)'#r'\n'
        t.lexer.lineno += len(t.value)
        self.lineno=t.lexer.lineno-1
        return t

    def t_error(self,t):
        print "Illegal character '%s'" % t.value[0]
        t.lexer.skip(1)

    t_LETTERS      = r'[a-zA-Z]+'
    t_DIGITS       = r'[0-9]+'
    t_SIGN        = r'[+-]'
    t_PERIOD      = r'\.'
    t_SEMICOLON   = r'\;'
    t_EQUALS      = r'\='
    t_OPENPAREN   = r'\('
    t_CLOSEPAREN  = r'\)'
    t_OPENBRACE   = r'\{'
    t_CLOSEBRACE  = r'\}'
    t_OPENBRACKET = r'\['
    t_CLOSEBRACKET= r'\]'
    t_COLON       = r'\:'
    t_BACKSLASH   = r'\\'
    t_QUOTE       = r'"'
    t_COMMA       = r','
    #t_LINEEND     = r'(\n|\r|\n\r|\r\n)'
    t_WHITESPACE  = r'[\t ]+'
    t_UNDERSCORE  = r'_'
    t_SLASH       = r'/'
    t_OTHER       = r'[\~\@\^\*\&\$>\'<\?\%]'

    def __init__(self,**kwargs):
        if 'lextab' not in kwargs:
            import getpass
            kwargs['lextab']="tmpcache_optimizedlextab_"+getpass.getuser()
        if 'optimize' not in kwargs:
            kwargs['optimize']=1
        self.lexer = lex.lex(module=self, **kwargs)
        self.lineno=0

    def test(self,data,log=None):
        self.lexer.lineno=1
        self.lineno=0
        self.lexer.input(data)
        while True:
            tok = self.lexer.token()
            if not tok: break
            if log:
                log.write(repr(tok)+'\n')
            else:
                print tok
        if log:
            log.close()

class CDLParser:
    precedence = (
        #('right', 'BACKSLASH'),
        #('left', 'TIMES', 'DIVIDE'),
    )
    def __init__(self,alexer=None,**kwargs):
        if alexer==None:
            alexer=CDLLexer()
        self.lexer=alexer
        # proper format checks
        self.lastVariableName=None
        self.variables=set()
        self.dimensions=set()
        self.globalname="___global"
        self.tokens=self.lexer.tokens
        self.parser=yacc.yacc(module=self,start='CDL',write_tables=0,debug=0,**kwargs)

    @property
    def lineno(self):
        return self.lexer.lineno

    def parse(self,data):
        return self.parser.parse(data)

    def p_error(self,p):
        print 'Error at ',p

    def p_CDL(self,p):
        'CDL : LETTERS WHITESPACE varname ignored OPENBRACE eol dimensions_declarations variables_declarations CLOSEBRACE eof'
        p[0]={'dimensions':p[7],
            'variables':p[8],
            'attributes':p[8][self.globalname]['attributes']}
        del p[0]['variables'][self.globalname]

    def p_eof(self,p):
        '''eof : ignored
               | eof eol ignored'''
        p[0]=None

    def p_variables_declarations(self,p):
        '''variables_declarations : ignored LETTERS COLON ignored eol
                                  | variables_declarations variable_declaration
                                  | variables_declarations variable_attribute_declaration
                                  | variables_declarations global_attribute_declaration
                                  | variables_declarations ignored eol'''
        if p[1] is None:
            if p[2]!='variables':
                raise
            p[0]=OrderedDict([(self.globalname,{'attributes':OrderedDict()} )])
        else:
            p[0]=p[1]
            if p[2] is None:
                pass
            elif len(p[2])==2:#variable declaration
                p[0][p[2][0]]=p[2][1]
            else:
                if p[2][1] in p[0][p[2][0]]['attributes']:
                    raise RuntimeError('Duplicate attribute %s for variable %s on line %i. Error' % (p[2][1],p[2][0],self.lineno))
                p[0][p[2][0]]['attributes'][p[2][1]]=p[2][2]


    def p_dimensions_declarations(self,p):
        '''dimensions_declarations : ignored LETTERS COLON ignored eol
                                   | dimensions_declarations dimension_declaration
                                   | dimensions_declarations ignored eol'''
        if p[1] is None:
            if p[2]!='dimensions':
                raise
            p[0]=OrderedDict()
        else:
            p[0]=p[1]
            if p[2]!=None:
                p[0].update(p[2])

    def p_dimension_declaration(self,p):
        '''dimension_declaration : ignored varname ignored EQUALS ignored DIGITS ignored SEMICOLON ignored eol
                                 | ignored varname ignored EQUALS ignored LETTERS ignored SEMICOLON ignored eol'''
        if p[2] in self.dimensions:
            raise RuntimeError('Dimension '+p[2]+' already in CDL. Duplicate on line %i' % self.lineno)
        self.dimensions.add(p[2])
        try:
            p[0]=OrderedDict([(p[2],int(p[6]))])
        except Exception,e:
            if p[6]!='UNLIMITED':
                raise
            p[0]=OrderedDict([(p[2],-1)])

    def p_empty(self,p):
        'empty :'
        p[0]=None

    def p_number(self,p):
        '''number : fixednumber
                  | floatingnumber'''
        p[0]=p[1]

    def p_fixednumber(self,p):
        '''fixednumber : intnumber
                       | intnumber LETTERS'''
        if len(p)==2:
            p[0]=numpy.array(int(p[1]),dtype='int32')
        elif len(p)==3:
            assert(p[2]=='L')
            p[0]=numpy.array(long(p[1]),dtype='int64')

    def p_floatingnumber(self,p):
        '''floatingnumber : LETTERS
                          | SIGN LETTERS
                          | intnumber PERIOD
                          | intnumber PERIOD exponent
                          | intnumber PERIOD DIGITS
                          | intnumber PERIOD DIGITS exponent
                          | intnumber PERIOD LETTERS
                          | intnumber PERIOD exponent LETTERS
                          | intnumber PERIOD DIGITS LETTERS
                          | intnumber PERIOD DIGITS exponent LETTERS'''
        t='double'
        endpt=len(p)
        if(p[endpt-1]=='f'): #float finite
            endpt=endpt-1
            t='float'
        x=''.join(p[1:endpt])
        if t!='float' and len(p) in (2,3) and x[-1]=='f': #float nonfinite
            x=x[:-1]
            t='float'
        p[0]=numpy.array(float(x),dtype=t)

    def p_exponent(self,p):
        'exponent : LETTERS intnumber'
        assert(p[1]=='e')
        p[0]=''.join(p[1:])

    def p_intnumber(self,p):
        '''intnumber : DIGITS
                     | SIGN DIGITS'''
        p[0]=''.join(p[1:])

    def p_commentcontent(self,p):
        '''commentcontent : stringcontent
                          | commentcontent QUOTE stringcontent'''
        p[0]=''.join(p[1:])

    def p_ignored(self,p):
        '''ignored : empty
                   | ignored WHITESPACE'''
        p[0]=None

    def p_eol(self,p):
        '''eol :  LINEEND
               |  SLASH SLASH commentcontent LINEEND'''
        p[0]=None

    def p_complete_dimslist(self,p):
        '''complete_dimslist : OPENPAREN ignored CLOSEPAREN
                             | OPENPAREN dimslist CLOSEPAREN'''
        if p[2] is None:
            p[0]=tuple()
        else:
            p[0]=tuple(p[2])

    def p_varname(self,p):
        '''varname : LETTERS
                   | UNDERSCORE
                   | varname UNDERSCORE
                   | varname LETTERS
                   | varname DIGITS'''
        p[0]=''.join(p[1:])

    def p_stringvalue(self,p):
        'stringvalue : QUOTE stringcontent QUOTE'
        p[0]=p[2]

    def p_stringcontent(self,p):
        '''stringcontent : empty
                         | stringcontent BACKSLASH QUOTE
                         | stringcontent BACKSLASH BACKSLASH
                         | stringcontent LETTERS
                         | stringcontent WHITESPACE
                         | stringcontent UNDERSCORE
                         | stringcontent PERIOD
                         | stringcontent COLON
                         | stringcontent SEMICOLON
                         | stringcontent EQUALS
                         | stringcontent OPENPAREN
                         | stringcontent CLOSEPAREN
                         | stringcontent OPENBRACE
                         | stringcontent CLOSEBRACE
                         | stringcontent OPENBRACKET
                         | stringcontent CLOSEBRACKET
                         | stringcontent COMMA
                         | stringcontent SIGN
                         | stringcontent SLASH
                         | stringcontent DIGITS
                         | stringcontent OTHER'''
        if p[1] is None:
            p[0]=''
        elif p[2]=='\\':
            p[0]=p[1]+p[3]
        else:
            p[0]=''.join(p[1:])

    def p_vectorvalue(self,p):
        'vectorvalue : vectorcontent'
        p[0]=p[1]

    def p_vectorcontent(self,p):
        '''vectorcontent : ignored number ignored
                         | vectorcontent COMMA ignored number ignored'''
        if p[1] is None:
            p[0]=p[2]
        else:
            p[0]=p[1]
            if len(p[0].shape)==0:
                p[0]=p[0].reshape([1])
            x=p[4].reshape([1])
            p[0]=numpy.concatenate([p[0],x],0)

    def p_typename(self,p):
        'typename : LETTERS'
        p[0]=p[1]

    def p_dimslist(self,p):
        '''dimslist : ignored varname ignored
                    | dimslist COMMA ignored varname ignored'''
        if p[1] is None:
            p[0]=[]
            p[0].append(p[2])
        else:
            p[0]=p[1]
            p[0].append(p[4])

    def p_variable_value(self,p):
        '''variable_value : ignored stringvalue ignored 
                          | vectorvalue'''
        if p[1] is None:
            p[0]=p[2]
        else:
            p[0]=p[1]

    def p_attribute_declaration(self,p):
        'attribute_declaration : COLON varname ignored EQUALS variable_value SEMICOLON ignored eol'
        p[0]=(p[2],p[5])

    def p_global_attribute_declaration(self,p):
        'global_attribute_declaration : ignored attribute_declaration'
        p[0]=(self.globalname,p[2][0],p[2][1])

    def p_variable_attribute_declaration(self,p):
        'variable_attribute_declaration : ignored varname attribute_declaration'
        if not p[2] in self.variables:
            raise RuntimeError('Variable Attribute for Nonexistant variable '+p[2]+' in CDL. Error on line %i' % self.lineno)
        if p[2] != self.lastVariableName:
            raise RuntimeError('Out of order attribute. Is for variable '+p[2]+', but under variable '+self.lastVariableName+' in CDL. Error on line %i' % self.lineno)            
        p[0]=(p[2],p[3][0],p[3][1])

    def p_variable_declaration(self,p):
        '''variable_declaration : ignored typename WHITESPACE varname ignored complete_dimslist ignored SEMICOLON ignored eol
                                | ignored typename WHITESPACE varname ignored SEMICOLON ignored eol'''
        if p[4] in self.variables:
            raise RuntimeError('Variable '+p[4]+' already in CDL. Duplicate on line %i' % self.lineno)
        self.variables.add(p[4])
        self.lastVariableName=p[4]
        p[0]=(p[4],{'type':p[2],'dimensions':p[6] if (len(p)>9) else tuple(),'attributes':OrderedDict()})

def CDL2dict(parm):
    p=CDLParser()
    if isinstance(parm,basestring):
        if os.path.exists(parm):
            parm=file(parm)
        else:
            return p.parse(parm)
    return p.parse(parm.read())

import functools
from collections import namedtuple
typecodes={
    'byte':'b',
    'char':'c',
    'short':'i2',
    'int':'i4',
    'long':'i8',
    'double':'d',
    'float':'f'
}

def getkey(d,key):
    return d[key]

def returnval(v):
    return v

def PseudoNetCDF4Variable(_d):
    d=_d.copy()
    attrs=d['attributes'].copy()
    if '_FillValue' in attrs:
        d['fill_value']=attrs['_FillValue']
        del attrs['_FillValue']
    elif 'type' in d:
        if d['type'] in ['char']:
            d['fill_value']=u'\x00'
        else:
            d['fill_value']=None
    d['ncattrs']=attrs.keys
    d['getncattr']=functools.partial(getkey,attrs)
    d.update(attrs)
    del d['attributes']
    remov=[]
    for f in d.keys():
        if f.startswith('_'):
            remov.append(f)
    for f in remov:
        del d[f]
    if 'type' in d:
        d['dtype']=typecodes[d['type']]
    return namedtuple('PseudoNetCDF4Variable',' '.join(d.keys()))(**d)

class PseudoNetCDF4Dimension(object):
    def __init__(self,length):
        self.length=length

    def isunlimited(self):
        return self.length is None or self.length<0

    def __len__(self):
        if self.isunlimited():
            return 0
        return self.length

def CDL2pseudonetCDF4(parm):
    d=CDL2dict(parm)
    for k,v in d['dimensions'].items():
        d['dimensions'][k]=PseudoNetCDF4Dimension(v)
    for k,v in d['variables'].items():
        d['variables'][k]=PseudoNetCDF4Variable(v)
    return PseudoNetCDF4Variable(d)

if __name__ == '__main__':
    m=CDLLexer()
    import json
    p=None#CDLParser()
    bindings=None
    nextBindings=False
    for x in range(1,len(sys.argv)):
        if sys.argv[x]=='-b':
            nextBindings=True
            continue
        if nextBindings:
            bindings=sys.argv[x].split(',')
            nextBindings=False
            continue
        f=file(sys.argv[x])
        m.test(f.read(),file('lexerout_'+(sys.argv[x].split('/')[-1]),'w'))
        f=file(sys.argv[x])
        if p is None:
            p=CDLParser()
        d=p.parse(f.read())
        if bindings!=None:
            b=[]
            for f in bindings:
                if f in d['variables']:
                    if 'dpl_py_binding' in d['variables'][f]['attributes']:
                        b.append( d['variables'][f]['attributes']['dpl_py_binding'])
                else:
                    print f,'not in template'
            print b
            bindings=None
        #json.dump(d,file((sys.argv[x].split('/')[-1])+'.json','w'),indent=4,separators=(',', ': '))
