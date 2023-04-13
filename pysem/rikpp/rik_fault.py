# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
RIK fault definition
"""
# Required modules
import numpy as np
from sem_fault import SEM3Dfault
from rik_pp_lib import classproperty, parseRIKrates
from scipy.integrate import cumtrapz


class RIK2DFault(SEM3Dfault):
    __nL = 1
    __nW = 1
    __nt = 0
    __dt = 0.0
    __HypoFile = ''
    __HypoAlongSDDepthKm = np.empty((3,), dtype=np.float64)
    __setHypoFile = False
    __setSlipFile = False
    __setMomentRateFile = False
    __setSlipRateFile = False
    __HypoFile = ''
    __SlipFile = ''
    __MomentRateFile = ''
    __SlipRateFile = ''
    __HypoDepthm = 0.0
    __SlipGridAlongS = np.empty((1, 1), dtype=np.float64)
    __SlipGridAlongD = np.empty((1, 1), dtype=np.float64)
    __SlipAlongSD = np.empty((1, 1), dtype=np.float64)
    __SlipGridAlongSD = np.empty((1, 1, 1), dtype=np.float64)
    __MomentGridAlongSD = np.empty((1, 1, 1), dtype=np.float64)
    __MomentRateGridAlongSD = np.empty((1, 1, 1), dtype=np.float64)
    __SlipRateGridAlongSD = np.empty((1, 1, 1), dtype=np.float64)

    def __init__(self, **kwargs):
        super(RIK2DFault, self).__init__()
        if 'nL' in kwargs.keys():
            RIK2DFault.__nL = kwargs['nL']
        if 'nW' in kwargs.keys():
            RIK2DFault.__nW = kwargs['nW']
        if 'nt' in kwargs.keys():
            RIK2DFault.__nt = kwargs['nt']

        if 'HypoAlongSDDepthKm' in kwargs.keys():
            RIK2DFault.__HypoAlongSDDepthKm = kwargs['HypoAlongSDDepthKm']

    @classproperty
    def nL(cls):
        return cls.__nL

    @classproperty
    def nW(cls):
        return cls.__nW

    @classproperty
    def nt(cls):
        return cls.__nt

    @classproperty
    def dt(cls):
        return cls.__dt

    @nL.setter
    def SetnL(cls, nL: int) -> None:
        cls.__nL = nL

    @nW.setter
    def SetnW(cls, nW: int) -> None:
        cls.__nW = nW

    @nt.setter
    def Setnt(cls, nt: int) -> None:
        cls.__nt = nt

    @dt.setter
    def Setnt(cls, dt: float) -> None:
        cls.__dt = dt

    @classproperty
    def setHypoFile(cls):
        return cls.__setHypoFile

    @setHypoFile.setter
    def SetsetHypoFile(cls, setHypoFile: bool):
        cls.__setHypoFile = setHypoFile

    @classproperty
    def HypoFile(cls):
        return cls.__HypoFile

    @HypoFile.setter
    def SetHypoFile(cls, HypoFile: str):
        cls.__HypoFile = HypoFile
        cls.setHypoFile = True

    @classproperty
    def HypoDepthKm(cls):
        return cls.__HypoAlongSDDepthKm[-1]

    @HypoDepthKm.setter
    def SetHypoDepthKm(cls, depth_km: float):
        cls.__HypoAlongSDDepthKm[-1] = depth_km
        print("Hypocenter location: depth: {} km\n\n".format(*
              tuple(cls.HypoAlongSDDepthKm[-1])))

    @classproperty
    def HypoDepthm(cls):
        return cls.__HypoDepthm

    @HypoDepthm.setter
    def SetHypoDepthm(cls, HypoDepthm: float):
        if HypoDepthm:
            cls.__HypoDepthm = HypoDepthm
        else:
            cls.__HypoDepthm = cls.HypoAlongSDDepthKm[-1]*1.0e3
        print("Hypocenter location: depth: {} m\n\n".format(cls.HypoDepthm))

    @classproperty
    def HypoAlongSD(cls):
        return cls.__HypoAlongSDDepthKm[:-1]

    @HypoAlongSD.setter
    def SetHypoAlongSD(cls, HypoDict: dict) -> None:
        HypoFile = HypoDict["HypoFile"]
        if not cls.setHypoFile:
            cls.HypoFile = HypoFile
        cls.HypoAlongSDDepthKm[0] = np.genfromtxt(cls.HypoFile, usecols=0)
        cls.HypoAlongSDDepthKm[1] = np.genfromtxt(cls.HypoFile, usecols=1)
        print("Hypocenter location: along-length: {} km - along-dip: {}\n\n".format(*
              tuple(cls.HypoAlongSDDepthKm[:-1])))

    @classproperty
    def HypoAlongSDDepthKm(cls):
        return cls.__HypoAlongSDDepthKm

    @HypoAlongSDDepthKm.setter
    def SetHypoAlongSDDepthKm(cls, HypoDict: dict) -> None:
        HypoFile = HypoDict["HypoFile"]
        HypoDepthKm = HypoDict["HypoDepthKm"]
        if not cls.setHypoFile:
            cls.HypoFile = HypoFile
        cls.HypoAlongSDDepthKm[0] = np.genfromtxt(cls.HypoFile, usecols=0)
        cls.HypoAlongSDDepthKm[1] = np.genfromtxt(cls.HypoFile, usecols=1)
        cls.HypoDepthKm = HypoDepthKm

        print("Hypocenter location: along-length: {} km - along-dip: {} km - depth: {} km\n\n".format(*
              tuple(cls.HypoAlongSDDepthKm)))

    @classproperty
    def setSlipFile(cls):
        return cls.__setSlipFile

    @setSlipFile.setter
    def SetsetSlipFile(cls, setSlipFile: bool):
        cls.__setSlipFile = setSlipFile

    @classproperty
    def SlipFile(cls):
        return cls.__SlipFile

    @SlipFile.setter
    def SetSlipFile(cls, SlipFile: str):
        cls.__SlipFile = SlipFile
        cls.setSlipFile = True

    @classproperty
    def SlipGridAlongS(cls):
        return cls.__SlipGridAlongS

    @classproperty
    def SlipGridAlongD(cls):
        return cls.__SlipGridAlongD

    @classproperty
    def GridAlongSD(cls):
        return cls.__SlipAlongSD

    @SlipGridAlongS.setter
    def SetSlipGridAlongS(cls, GridAlongS: np.float64) -> None:
        cls.__SlipGridAlongS = GridAlongS

    @SlipGridAlongD.setter
    def SetSlipGridAlongD(cls, GridAlongD: np.float64) -> None:
        cls.__SlipGridAlongD = GridAlongD

    @GridAlongSD.setter
    def SetGridAlongSD(cls, SlipFileDict: dict) -> None:
        SlipFile = SlipFileDict['SlipFile']
        if not cls.setSlipFile:
            cls.SlipFile = SlipFile
        cls.SlipGridAlongS = np.genfromtxt(cls.SlipFile,
                                           usecols=1).reshape(cls.nL,
                                                              cls.nW,
                                                              order='F')
        cls.SlipGridAlongD = np.genfromtxt(cls.SlipFile,
                                           usecols=0).reshape(cls.nL,
                                                              cls.nW,
                                                              order='F')

    @classproperty
    def setSlipRateFile(cls):
        return cls.__setSlipRateFile

    @setSlipRateFile.setter
    def SetsetSlipRateFile(cls, setSlipRateFile: bool):
        cls.__setSlipRateFile = setSlipRateFile

    @classproperty
    def setMomentRateFile(cls):
        return cls.__setMomentRateFile

    @setMomentRateFile.setter
    def SetsetMomentRateFile(cls, setMomentRateFile: bool):
        cls.__setMomentRateFile = setMomentRateFile

    @classproperty
    def MomentRateFile(cls):
        return cls.__MomentRateFile

    @classproperty
    def SlipRateFile(cls):
        return cls.__SlipRateFile

    @MomentRateFile.setter
    def SetMomentRateFile(cls, MomentRateFile: str) -> None:
        cls.__MomentRateFile = MomentRateFile
        cls.setMomentRateFile = True

    @SlipRateFile.setter
    def SetSlipRateFile(cls, SlipRateFile: str) -> None:
        cls.__SlipRateFile = SlipRateFile
        cls.setSlipRateFile = True

    @classproperty
    def MomentRateGridAlongSD(cls):
        return cls.__MomentRateGridAlongSD

    @MomentRateGridAlongSD.setter
    def SetMomentRateGridAlongSD(cls, MomentRateFileDict: dict) -> None:
        MomentRateFile = MomentRateFileDict['MomentRateFile']
        if not cls.setMomentRateFile:
            cls.MomentRateFile = MomentRateFile
        cls.__MomentRateGridAlongSD = parseRIKrates(cls.MomentRateFile,
                                                    (cls.nL,
                                                     cls.nW,
                                                     cls.nt)
                                                    )
        cls.MomentGridAlongSD = cumtrapz(cls.MomentRateGridAlongSD,
                                         dx=cls.dt,
                                         initial=0.,
                                         axis=-1)

    @classproperty
    def SlipRateGridAlongSD(cls):
        return cls.__SlipRateGridAlongSD

    @SlipRateGridAlongSD.setter
    def SetSlipRateGridAlongSD(cls, SlipRateFileDict: dict) -> None:
        SlipRateFile = SlipRateFileDict['SlipRateFile']
        if not cls.setMomentRateFile:
            cls.SlipRateFile = SlipRateFile
        cls.__SlipRateGridAlongSD = parseRIKrates(cls.SlipRateFile,
                                                  (cls.nL,
                                                   cls.nW,
                                                   cls.nt)
                                                  )

        cls.SlipGridAlongSD = cumtrapz(cls.SlipRateGridAlongSD,
                                       dx=cls.dt,
                                       initial=0.,
                                       axis=-1)

    @classproperty
    def SlipGridAlongSD(cls):
        return cls.__SlipGridAlongSD

    @SlipGridAlongSD.setter
    def SetSlipGridAlongSD(cls, SlipGridAlongSD: np.float64) -> None:
        cls.__SlipGridAlongSD = SlipGridAlongSD

    @classproperty
    def MomentGridAlongSD(cls):
        return cls.__MomentGridAlongSD

    @MomentGridAlongSD.setter
    def SetMomentGridAlongSD(cls, MomentGridAlongSD: np.float64) -> None:
        cls.__MomentGridAlongSD = MomentGridAlongSD