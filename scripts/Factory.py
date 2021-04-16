from assign_clades import alignArgs, tmpNode
from modify_tree_according_to_exposure import TraverseImpl,ModifyImpl,Collect_Exposures_IMPL,SAMP_TRAIT_IMPL,Update_latlongs_IMPL

class Factory(object):
    __alignArgs = None
    __tmpNode = None
    __TraverseImp,__Collect_Exposures_IMPL,__SAMP_TRAIT_IMPL = None,None,None
    __ModifyImpl,__Update_latlongs_IMPL = None,None

    def __init__(self):
        pass

    def getAlignArgs(self):
        if Factory.__alignArgs is None:
            return alignArgs()
        else:
            return Factory.__alignArgs

    def getTmpNode(self):
        if Factory.__tmpNode is None:
            return tmpNode()
        else:
            return Factory.__tmpNode

    def getTraverse(self):
        if Factory.__TraverseImp is None:
            return TraverseImpl()
        else:
            return Factory.__TraverseImp

    def getcollectexposure(self):
        if Factory.__Collect_Exposures is None:
            return Collect_Exposures_IMPL()
        else:
            return Factory.__Collect_Exposures_IMPL

    def getsamptrait(self):
        if Factory.__SAMP_TRAIT_IMPL is None:
            return SAMP_TRAIT_IMPL()
        else:
            return Factory.__SAMP_TRAIT_IMPL

    def getmodify(self):
        if Factory.__ModifyImpl is None:
            return ModifyImpl()
        else:
            return Factory.__ModifyImpl

    def getupdatedlatlong(self):
        if Factory.__Update_latlongs_IMPL is None:
            return Update_latlongs_IMPL()
        else:
            return Update_latlongs_IMPL







Factory.getAlignArgs = staticmethod(Factory.getAlignArgs)
Factory.getTmpNode = staticmethod(Factory.getTmpNode)
Factory.getTraverse = staticmethod(Factory.getTraverse)
Factory.getcollectexposure = staticmethod(Factory.getcollectexposure)
Factory.getsamptrait = staticmethod(Factory.getsamptrait)
Factory.getmodify = staticmethod(Factory.getmodify)
Factory.getupdatedlatlong = staticmethod(Factory.getupdatedlatlong)
