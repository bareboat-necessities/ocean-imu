"""Positive finite range of the pinned powf at the two tuner exponents.

The proof follows newlib's actual log2/high-low/exp2 graph. It establishes
a range, not approximation accuracy or a correctly-rounded fractional root.
Exact rational bounds include normal/subnormal RNE error and contractions.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import shlex
import tempfile

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import target_qaxis_exp as Q
from tools.stability.ou3_alt_contraction import target_wpe_libm as W
from tools.stability.ou3_alt_contraction import target_scalar_profile as S
from tools.stability.ou3_alt_contraction import target_math_link as LINK

QUALIFICATION='OU3_ALT_PINNED_POWF_TUNER_FINITE_RANGE_V1'
SOURCE_SHA='ca8a3229cd7f4cf68012f215b25a287627b6ed92120cb2a65679c096536f8626'
MEMBER_SHA='9dff07e2ab4b86832b743f706b25431ec1df1ec3534296d367cf9088b746a2cb'
BASE_DOMAIN=(F(1,2**44),F(2**17))
EXPONENTS=(F(14380471,16777216),F(9586981,134217728))
OUTPUT_RANGE=(F(1,2**45),F(2**18))
U=F(1,2**24); ETA=F(1,2**150); DROP=F(1,2**11)
LWORDS=(0x3f19999a,0x3edb6db7,0x3eaaaaab,0x3e8ba305,0x3e6c3255,0x3e53f142)
L=tuple(map(W._word,LWORDS))
CP,CP_H,CP_L=map(W._word,(0x3f76384f,0x3f763800,0x369dc3a0))
DP_H,DP_L=map(W._word,(0x3f15c000,0x35d1cfdc))
LG2,LG2_H,LG2_L=map(W._word,(0x3f317218,0x3f317200,0x35bfbe8c))


def up(m): return F(m)*(1+U)+ETA
def add(a,b): return up(F(a)+F(b))
def mul(a,b): return up(F(a)*F(b))
def drop_error(m): return DROP*F(m)+4096*ETA
def round_error(m): return U*F(m)+ETA


def certificate():
    # Three exact mantissa branches from the source's bit selectors. Each
    # endpoint denotes an entire contiguous mantissa interval, not a sample.
    branches=((0x3f800000,0x3f9cc471,F(1),0),
              (0x3f9cc472,0x3fddb3d6,F(3,2),1),
              (0x3f6db3d7,0x3f7fffff,F(1),0))
    ratio=F(0)
    for lo,hi,bp,k in branches:
        a,b=W._word(lo),W._word(hi)
        assert bp/2<=a<=b<=2*bp # ax-bp exact by Sterbenz
        assert a+bp>F(9,5) and b+bp<F(13,4)
        ratio=max(ratio,abs((a-bp)/(a+bp)),abs((b-bp)/(b+bp)))
        # This bit expression is monotone within each fixed-exponent branch.
        th=[W._word(((word>>1)|0x20000000)+0x40000+(k<<21)) for word in (lo,hi)]
        assert 1<=min(th)<=max(th)<=4
    # u is exact, v=RN(1/RN(ax+bp)), s=RN(u*v).
    ratio_error=ratio*((1+U)**2/(1-U)-1)+16*ETA
    sm=ratio+ratio_error
    assert sm<F(1,9)
    sm=F(1,9); vm=up(1/(F(9,5)*(1-U)))
    # t_l=RN(ax-RN(t_h-bp)); preserve t_h+t_l=ax+bp+delta.
    tl_error=round_error(3)+round_error(5+round_error(3))
    residual=F(13,4)*(ratio_error+drop_error(sm))+sm*tl_error
    # u-sh*t_h-sh*t_l: charge every product and subtraction, including a
    # larger separate-rounding allowance for either fused MSUB choice.
    residual+=round_error(sm*4)+round_error(F(7,25)+mul(sm,4))
    residual+=round_error(sm*5)+round_error(add(F(7,25),mul(sm,4))+mul(sm,5))
    sl=mul(vm,residual)
    assert sl<F(1,4096)
    sl=F(1,4096)
    s2=mul(sm,sm); poly=abs(L[-1])
    for coefficient in reversed(L[:-1]): poly=add(abs(coefficient),mul(s2,poly))
    rm=add(mul(mul(s2,s2),poly),mul(sl,add(sm,sm)))
    th=add(add(3,s2),rm)
    # t_l=r-((trunc(RN(3+s_h²+r))-3)-s_h²): retain compensation.
    tlo=drop_error(th)+round_error(3+s2)+round_error(add(3,s2)+rm)
    tlo+=round_error(th+3)+round_error(th+3+s2)+round_error(rm+th+3+s2)
    assert tlo<F(1,512)
    tlo=F(1,512)
    um=mul(sm,th); vmag=add(mul(sl,th),mul(tlo,sm))
    ph=add(um,vmag) # truncation cannot increase magnitude
    pl=drop_error(ph)+round_error(um+vmag)+round_error(ph+um)
    pl+=round_error(vmag+ph+um)
    assert pl<F(1,2048)
    pl=F(1,2048)
    zh=mul(CP_H,ph); zl=add(add(mul(CP_L,ph),mul(pl,CP)),DP_L)
    # Input exponent n is [-44,17]. The third mantissa branch increments n;
    # at the upper endpoint 2^17 only j=0 is legal, so n never reaches 18.
    t1lo=-up(up(up(zh+zl))+44)
    t1hi=add(add(add(zh,zl),DP_H),17)
    assert t1lo>-45 and t1hi<18
    # t1 is the chopped three-add sum; bound t2 without discarding its
    # cancellation against the SAME n, dp_h, z_h and z_l.
    t2=drop_error(45)+round_error(zh+zl)+round_error(add(zh,zl)+DP_H)
    t2+=round_error(add(add(zh,zl),DP_H)+44)
    t2+=round_error(45+44)+round_error(45+44+DP_H)
    t2+=round_error(45+44+DP_H+zh)+round_error(zl+45+44+DP_H+zh)
    assert t2<F(1,32)
    t2=F(1,32)
    # Both actual binary32 exponent constants are positive, below 7/8,
    # distinct from powf's special 0,1,2,1/2 and huge-y branches.
    assert all(F(1,16)<y<F(7,8) and B.is_binary32(y) for y in EXPONENTS)
    ymax=F(7,8)
    pm=add(mul(drop_error(ymax),45),mul(ymax,t2))
    product_error=round_error(drop_error(ymax)*45)+round_error(ymax*t2)
    product_error+=round_error(pm)+round_error(ymax*45)+round_error(pm+ymax*45)
    zlo=-ymax*(45+t2)-product_error
    zhi=ymax*(18+t2)+product_error
    assert zlo>-40 and zhi<16
    # The bit add/mask in ef_pow.c gives sign(z)*floor(abs(z)+1/2).
    # Its integer n is therefore [-40,16], and |z-n|<=1/2. Subtracting
    # float(n) from p_h and readding p_l adds only these round errors.
    red=F(1,2)+round_error(ymax*45+40)+round_error(pm+ymax*45)+round_error(1)
    assert red<F(501,1000)
    red=F(501,1000)
    # t=trunc(RN(p_h+p_l)); hence p_l-(t-p_h) is a small residual.
    low=drop_error(up(red))+round_error(red)+round_error(red+(red+pm))
    low+=round_error(pm+red+(red+pm))
    u2=mul(red,LG2_H); v2=add(mul(low,LG2),mul(red,LG2_L))
    z2=add(u2,v2)
    assert z2<F(7,20)
    z2=F(7,20)
    # w=v-(z-u) uses the SAME computed z. FMA instead of RN(u)+v
    # removes one rounding and is covered by the extra u product error.
    wm=round_error(u2+v2)+round_error(z2+u2)+round_error(v2+z2+u2)
    wm+=round_error(u2)
    assert wm<F(1,2**20)
    wm=F(1,2**20)
    zz=mul(z2,z2); poly=abs(Q.P[-1])
    for coefficient in reversed(Q.P[:-1]): poly=add(abs(coefficient),mul(zz,poly))
    t1m=add(z2,mul(zz,poly))
    denominator=(2-t1m)*(1-U)-ETA
    assert denominator>F(3,2)
    remainder=add(up(mul(z2,t1m)/denominator),add(wm,mul(z2,wm)))
    radius=add(remainder,z2)
    core_lo=(1-radius)*(1-U)-ETA; core_hi=up(1+radius)
    assert F(1,2)<core_lo<core_hi<F(3,2)
    # The actual integer exponent-field addition is exact scaling: resulting
    # exponents [-41,16] remain normal. The scalbnf branch cannot be reached.
    actual_range=(F(1,2**41),F(3*2**15))
    assert OUTPUT_RANGE[0]<actual_range[0]<actual_range[1]<OUTPUT_RANGE[1]
    return {'qualification':QUALIFICATION,'base_domain':BASE_DOMAIN,
        'actual_exponents':EXPONENTS,'output_range':OUTPUT_RANGE,
        'proved_stronger_output_range':actual_range,
        'mantissa_ratio_bound':sm,'log2_high_interval':(F(-45),F(18)),
        'log2_low_absolute_bound':t2,'power_exponent_interval':(F(-40),F(16)),
        'exp2_reduced_argument_absolute_bound':z2,
        'exp2_reduced_output_interval':(core_lo,core_hi),
        'all_intermediates_finite':True,'return_strictly_positive':True,
        'exception_and_subnormal_output_branches_unreachable':True,
        'approximation_accuracy_qualified':False}


def audit(toolchain:Path,sdk:Path):
    compiler=toolchain/'bin/xtensa-esp32s3-elf-g++'
    ar=toolchain/'bin/xtensa-esp32s3-elf-ar'; od=toolchain/'bin/xtensa-esp32s3-elf-objdump'
    libm=Path(S.run([compiler,'-print-file-name=libm.a']).strip())
    if S.sha(libm)!=Q.LIBM_SHA256: raise ValueError('unqualified powf libm archive')
    import subprocess
    raw=subprocess.check_output([str(ar),'p',str(libm),'libm_a-ef_pow.o'])
    if hashlib.sha256(raw).hexdigest()!=MEMBER_SHA: raise ValueError('powf kernel bytes changed')
    source='#include <cmath>\nextern "C" float alt_pow_probe(float x,float y){return std::pow(x,y);}\n'
    flags=shlex.split((sdk/'flags/cpp_flags').read_text())+['-Os','-funroll-loops','-fno-finite-math-only']
    for ld in (sdk/'ld').glob('*.ld'):
        if re.search(r'\b(?:powf|__ieee754_powf)\s*=',ld.read_text()):
            raise ValueError('SDK powf symbol override requires qualification')
    with tempfile.TemporaryDirectory(prefix='ou3-pow-link-') as directory:
        tmp=Path(directory); cpp=tmp/'probe.cpp'; obj=tmp/'probe.o'; elf=tmp/'probe.elf'
        layout=tmp/'layout.ld'; mapping=tmp/'probe.map'; member=tmp/'pow.o'
        cpp.write_text(source); layout.write_text(LINK.LAYOUT); member.write_bytes(raw)
        S.run([compiler,*flags,'-ffile-prefix-map='+str(tmp)+'=/ou3-pow-link','-c',cpp,'-o',obj])
        linked=LINK._run([compiler,'-nostartfiles','-Wl,-e,alt_pow_probe','-Wl,--gc-sections',
            '-Wl,-T,'+str(layout),'-Wl,-T,'+str(sdk/'ld/esp32s3.rom.libgcc.ld'),
            '-Wl,-Map,'+str(mapping),obj,'-o',elf,'-lm','-lc','-lgcc'])
        for name in ('libm_a-ef_pow.o','libm_a-wf_pow.o'):
            if '('+name+')' not in mapping.read_text(): raise ValueError('powf namespace supplier changed')
        dis=S.run([od,'-dr',member]); linked_dis=S.run([od,'-d',elf])
        # Audit the actual normalized-log and reduced-exp instructions; the
        # source is built with contraction and the bound allows those fusions.
        operations={op:len(re.findall(r'\s'+re.escape(op)+r'\s',dis))
                    for op in ('add.s','sub.s','mul.s','madd.s','msub.s')}
        literals=S.run([od,'-s','-j','.literal.__ieee754_powf',member])
        for word in (*LWORDS,*Q.COEFFICIENT_WORDS,0x3f76384f,0x3f763800,0x369dc3a0):
            if word.to_bytes(4,'little').hex() not in literals:
                raise ValueError('powf coefficient literal missing')
        if not re.search(r'40002274\s+A\s+__divsf3',S.run([toolchain/'bin/xtensa-esp32s3-elf-nm',elf])):
            raise ValueError('powf division does not use audited ROM entry')
        report={'qualification':QUALIFICATION,'source_commit':Q.NEWLIB_COMMIT,
            'source_sha256':SOURCE_SHA,'libm_sha256':S.sha(libm),
            'powf_kernel_member_sha256':MEMBER_SHA,'kernel_arithmetic_instruction_counts':operations,
            'compiler_driver_sha256':S.sha(compiler.resolve()),'compiler_flags':flags,
            'sdk_cpp_flags_sha256':S.sha(sdk/'flags/cpp_flags'),
            'standalone_powf_namespace_resolution_qualified':True,
            'powf_division_resolves_to_audited_ROM_entry':True,
            'standalone_runtime_complete':False,
            'non_math_runtime_linker_diagnostic':linked.stderr.replace(str(toolchain),'/toolchain').strip(),
            'whole_firmware_link_and_callgraph_qualified':False,
            'rounding_premise':'FCR.RM=0 throughout the admitted execution'}
    report['range_certificate']=certificate()
    return report


def profile_correspondence():
    path=Path(__file__).resolve().parents[1]/'ou3_alt_target_powf_range.json'
    report=json.loads(path.read_text())
    target=json.loads((path.parent/'ou3_alt_target_arithmetic_audit.json').read_text())
    if (report.get('qualification')!=QUALIFICATION or report.get('libm_sha256')!=Q.LIBM_SHA256
            or report.get('powf_kernel_member_sha256')!=MEMBER_SHA
            or any(report.get(k)!=target.get(k) for k in
                   ('compiler_driver_sha256','compiler_flags','sdk_cpp_flags_sha256'))):
        raise ValueError('stale pinned powf qualification')
    certificate()
    return {'pinned_tuner_powf_positive_finite_range_closed':(
        report.get('standalone_powf_namespace_resolution_qualified') is True
        and W.profile_correspondence()['pinned_scalar_profile_attached']),
        'base_domain':BASE_DOMAIN,'actual_exponents':EXPONENTS,'output_range':OUTPUT_RANGE,
        'whole_firmware_compiler_and_link_correspondence_closed':False}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('toolchain','sdk','output'): p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args(); result=audit(a.toolchain,a.sdk)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True,default=str)+'\n')
    print(json.dumps(profile_correspondence(),default=str))


if __name__=='__main__': main()
