from native_trimers import *
report=json.loads((R/'results/paralog_profile_summary.json').read_text());errors=[]
for lam in [1.,4.]:
 row=next(x for x in report if x['lam']==lam and x['strain']=='WT');v=np.log(row['parameters'])
 a,qa=simulate('WT',lam,v,method='LSODA',rtol=2e-8);b,qb=simulate('WT',lam,v,method='BDF',rtol=2e-8)
 err=float(np.max(np.abs(a[:,4:]-b[:,4:])));assert err<2e-5
 errors.append(dict(lam=lam,max_BDF_LSODA_difference=err,LSODA=qa,BDF=qb))
(R/'results/native_independent_verification.json').write_text(json.dumps(errors,indent=2));print(json.dumps(errors,indent=2))
