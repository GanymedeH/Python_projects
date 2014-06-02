import matplotlib.pyplot as mpl
import pickle
import time
import os
import numpy as np





def qpick(fname, mode="r", data=None):
    import cPickle
    f=open(fname,mode)
    if mode == "r":
        data = cPickle.load(f)
    elif mode == "w":
        cPickle.dump(data, f)
    f.close()
    return data









def plot_vdata(filename="recvs.pickle", step=1, plottype="theta"):
    runID = int(time.time())
    vf = open(filename,"r")
    vvals = pickle.load(vf)
    vf.close()

    # pick these by trial and error
    pt1 = 1
    pt2 = 4
    pt3 = 7
    # plots v value vs arc ID, both sorted and unsorted (sorted destroys ID# across time steps). tracks 3 pts above
    itern = 0
    for vkt in vvals:
        if itern % step != 0:
            itern += 1
            continue
        if itern < 10:
            itern = "000"+str(itern)
        elif itern < 100:
            itern = "00"+str(itern)
        elif itern < 1000:
            itern = "0"+str(itern)
        cnt = 1
        for vt in vkt:
            if cnt > 5:
                break
            fcnt = 1
            for vf in vt:
                mpl.figure(fcnt)
                mpl.subplot(2,3,cnt)
                mpl.plot(vf, color='k', hold=True)
                mpl.plot([pt1],[vf[pt1]],marker='^',mfc='c',mec='y',ms=7,mew=1)
                mpl.plot([pt2],[vf[pt2]],marker='o',mfc='r',mec='y',ms=7,mew=1)
                mpl.plot([pt3],[vf[pt3]],marker='s',mfc='y',mec='g',ms=7,mew=1)
                mpl.xticks(size="small")
                if plottype == "theta":
                    mpl.yticks(size="x-small")
                    if cnt == 2:
                        mpl.title("Quadrant "+str(fcnt)+" value approximation coefficient vs arc ID")
                    mpl.xlabel("time "+str(cnt))
                elif plottype == "del_theta":
##                    mpl.ylim(-15, 15)
                    mpl.yticks(size="x-small")
                    if cnt == 2:
                        mpl.title("Quadrant "+str(fcnt)+\
                                  " value coefficient delta over iterations vs arc ID @ iteration "+itern)
                    mpl.xlabel("time "+str(cnt))
                    # don't make a sorted plot (wouldn't make sense)
                    fcnt += 1
                    continue

                evf = list(enumerate(vf))
                vfs = sorted(evf, key=lambda x: x[1])
                idx1 = vfs.index([(x,y) for x,y in vfs if x == pt1][0])
                idx2 = vfs.index([(x,y) for x,y in vfs if x == pt2][0])
                idx3 = vfs.index([(x,y) for x,y in vfs if x == pt3][0])

                mpl.figure(len(vt)+fcnt+1)
                mpl.subplot(2,3,cnt)
                mpl.plot([y for x,y in vfs], color='k', hold=True)
                mpl.plot([idx1],[vfs[idx1][1]],marker='^',mfc='c',mec='y',ms=7,mew=1)
                mpl.plot([idx2],[vfs[idx2][1]],marker='o',mfc='r',mec='y',ms=7,mew=1)
                mpl.plot([idx3],[vfs[idx3][1]],marker='s',mfc='y',mec='g',ms=7,mew=1)
                mpl.xticks(size="small")
                ytickloc, yticklabel = mpl.yticks()
                mpl.yticks(size="x-small")
                if cnt == 2:
                    mpl.title("Sorted quadrant"+str(fcnt)+\
                              " value approximation coefficient @ iteration "+itern+"; arc ID not preserved", \
                              size="medium")
                mpl.xlabel("time "+str(cnt))
                fcnt += 1
            cnt += 1

        for ff in xrange(1,len(vt)+1):
            mpl.figure(ff)
            if plottype == "theta":
                if not os.path.exists("sorted_theta"):
                    os.mkdir("sorted_theta")
                if not os.path.exists("unsorted_theta"):
                    os.mkdir("unsorted_theta")
                mpl.savefig("./unsorted_theta/ADP_run_"+str(runID)+"_quad"+str(ff)+"_iter"+itern+\
                            "_v_vs_ID_unsorted.png",format="png")
##                mpl.clf()
                mpl.figure(len(vt)+ff+1)
                mpl.savefig("./sorted_theta/ADP_run_"+str(runID)+"_quad"+str(ff)+"_iter"+itern+\
                            "_v_vs_ID_sorted.png",format="png")
##                mpl.clf()
            elif plottype == "del_theta":
                if not os.path.exists("del_theta"):
                    os.mkdir("del_theta")
                mpl.savefig("./del_theta/ADP_run_"+str(runID)+"_quad"+str(ff)+"_iter"+itern+\
                            "_delv_vs_ID_unsorted.png",format="png")
##                mpl.clf()
        mpl.close("all")
        itern = int(itern)+1
    return






def plot_time_benchmarks(filename="cplex_solve_time_trials.txt"):
    runID = int(time.time())
    bf = open(filename,"r")
    loopcnt = 0
    cumdmts = 0
    firsttime=True
    mts = []
    delmts = []
    tott = [0]
    for bl in bf:
        if "main loop" in bl:
            for tok in bl.split():
                try:
                    mts.append(float(tok))
                    if firsttime:
                        firsttime=False
                    else:
                        delmts.append(float(tok)-oldmts)
                        cumdmts += delmts[-1]
                    oldmts = float(tok)
                    tott.append(tott[-1] + mts[-1])
                    loopcnt += 1
                except ValueError:
                    pass
    bf.close()
    mpl.figure(1)
    mpl.clf()
    mpl.plot(mts)
    mpl.title("Runtime vs. iteration #; Total runtime "+str(tott[-1])+" sec")
    mpl.ylabel("sec")
    mpl.xlabel("iteration #")
    mpl.savefig("ADP_run_"+str(runID)+"_mainloop_runtime_vs_iters.png",format="png")
    mpl.figure(2)
    mpl.clf()
    mpl.plot(tott)
    mpl.title("Cumulative runtime; Avg iteration length "+str(tott[-1]/float(loopcnt))+" sec")
    mpl.ylabel("sec")
    mpl.xlabel("iteration #")
    mpl.savefig("ADP_run_"+str(runID)+"_cumulative_runtime.png",format="png")
    mpl.clf()
    mpl.plot(delmts)
    mpl.title("Delta runtime vs. iteration #; Avg delta "+str(cumdmts/float(loopcnt))+" sec")
    mpl.ylabel("sec")
    mpl.xlabel("iteration #")
    mpl.savefig("ADP_run_"+str(runID)+"_delta_runtime.png",format="png")
    return





def plot_saa_benchmarks(contrupfile="an200_013_vcontrol.pickle", contrdownfile=None, \
                        testfiles=["an200_013_vtest.pickle"], legnd=None, scalefctr="2", \
                        combpng="an200_013_ip_vs_lp_saa_benchmarks_200mod2.png", seppngs=None, \
                        detailpngson=False, polouts=["an200_013_policyInstallPicks.pickle"]):
    vcu=qpick(contrupfile)
    if contrdownfile != None:
        vcd=qpick(contrdownfile)
    vt = []
    for tfn in testfiles:
        vt.append(qpick(tfn))
        vt[-1] = [i for i in vt[-1] if i[0] != None]

    if polouts != None:
        for k in xrange(len(polouts)):
            polpick=qpick(polouts[k])
            polpick=[i for i in polpick if i[0] != None]
            for i in xrange(len(polpick)-1):
                for j in xrange(len(polpick[i])):
                    if polpick[i+1][j] != polpick[i][j]:
                        print "policy change: test "+str(k)+" forecast "+str(j)+\
                              " from adp iteration "+str(i)+" to "+str(i+1)
            print "\n"
    for k in xrange(len(vt)):
        for i in xrange(len(vt[k])-1):
            for j in xrange(len(vt[k][i])):
                if vt[k][i+1][j] != vt[k][i][j]:
                    print "LPSAA change: test "+str(k)+" forecast "+str(j)+\
                          " from adp iteration "+str(i)+" to "+str(i+1)
        print "\n"

##    mpl.plot(range(len(vt)), np.ones(len(vt))*vcu, 'k', linewidth=2)
    mpl.figure()
    mpl.plot(range(len(vt[0])), np.ones(len(vt[0]))*np.mean(vcu), 'k', linewidth=2)
    mpl.hold("on")
    if contrdownfile != None:
        mpl.plot(range(len(vt[0])), np.ones(len(vt[0]))*np.mean(vcd), 'k', linewidth=2)
    for j in xrange(len(vt)):
        mpl.plot(range(len(vt[j])), [np.mean(i) for i in vt[j]], linewidth=2)
    mpl.xlabel("ADP training iteration count / "+scalefctr)
    mpl.ylabel("SAA objective value")
    mpl.title("SAA benchmark progression")
    if legnd == None:
        legnd = ["mu(LPSAA"+str(i)+")" for i in xrange(len(vt))]
        legnd.insert(0, "mu(IPSAA)-t0")
        if contrdownfile != None:
            legnd.insert(1, "mu(IPSAA)-tT")

    mpl.legend(legnd,loc=0)
    mpl.savefig(combpng,format="png")

    if seppngs == None:
        seppngs = ["an200_013_lp-ip_saa_diffs_200mod2_"+str(i)+".png" for i in xrange(len(vt))]

    optpct = []
    for k in xrange(len(vt)):
        mpl.figure()

        optgap = []
        optpct.append([])
        for i in xrange(len(vt[k])):
            optgap.append([vt[k][i][j] - vcd[j] for j in xrange(len(vcd))])
            optpct[-1].append([optgap[i][j] / vcd[j] for j in xrange(len(vcd))])
        ogmean = [np.mean(j) for j in optpct[-1]]

        mpl.plot(range(len(vt[k])), ogmean, 'k', linewidth=2)
##        mpl.hold("on")
##        mpl.plot(range(len(vt)), [j[4] for j in optpct], 'g', linewidth=2)
##        mpl.plot(range(len(vt)), [j[8] for j in optpct], 'r', linewidth=2)
##        mpl.plot(range(len(vt)), [j[12] for j in optpct], 'c', linewidth=2)

        mpl.xlabel("ADP training iteration count / "+scalefctr)
        mpl.ylabel("(LPSAA obj've - IPSAA(tT) obj've) / IPSAA obj've(tT)")
        mpl.title("Progression of LP/IP SAA optimality gap - percentages of IP SAA at time T")
##        mpl.legend(["mu(LP-IP)","LP-IP #4","LP-IP #8","LP-IP #12"],loc=5)
        mpl.savefig(seppngs[k],format="png")

    if detailpngson:
        for i in xrange(len(vt)):
            mpl.figure()

            mpl.subplot(2,3,1)
            mpl.plot(range(len(vt[i])), [j[3] for j in optpct[i]], 'k', linewidth=2)
            mpl.tick_params(labelsize=10)
            mpl.ylabel("Optimality gap (%)")

            mpl.subplot(2,3,2)
            mpl.plot(range(len(vt[i])), [j[6] for j in optpct[i]], 'k', linewidth=2)
            mpl.tick_params(labelsize=9)
            mpl.title("Selected optimality gaps\n\n",fontsize=16)

            mpl.subplot(2,3,3)
            mpl.plot(range(len(vt[i])), [j[12] for j in optpct[i]], 'k', linewidth=2)
            mpl.tick_params(labelsize=10)

            mpl.subplot(2,3,4)
            mpl.plot(range(len(vt[i])), [j[13] for j in optpct[i]], 'k', linewidth=2)
            mpl.tick_params(labelsize=10)
            mpl.ylabel("Optimality gap (%)")

            mpl.subplot(2,3,5)
            mpl.plot(range(len(vt[i])), [j[14] for j in optpct[i]], 'k', linewidth=2)
            mpl.tick_params(labelsize=9)
            mpl.xlabel("ADP training iteration count / "+scalefctr)

            mpl.subplot(2,3,6)
            mpl.plot(range(len(vt[i])), [j[19] for j in optpct[i]], 'k', linewidth=2)
            mpl.tick_params(labelsize=9)

            mpl.savefig("selectedsfrom_"+seppngs[i],format="png")
    return





if __name__ == "__main__":
##    plot_vdata(filename="recvs.pickle", step=4, plottype="theta")
##    print "v's plotted"
##    plot_vdata(filename="delvs.pickle", step=4, plottype="del_theta")
##    print "delv's plotted"
##    plot_time_benchmarks(filename="cplex_solve_time_trials.txt")
##    print "benchmarks plotted"

    tfiles = ["an200_013_vtest.pickle"]
    # even if None, it will still insert a generic legend
    leg = ["IP@time0", "IP@timeT","LP"]
    # even if None, it will still plot separate graphs for each test confiiguration
    separates = ["an200_013_test_mu_300mod5_0p03eps_0p03mcclain_ss.png"]
    pfiles = ["an200_013_policyInstallPicks.pickle"]
    plot_saa_benchmarks(contrupfile="an200_013_vcontrol_T0.pickle", \
                        contrdownfile="an200_013_vcontrol_Tfinal.pickle", testfiles=tfiles, legnd=leg, \
                        scalefctr="5", combpng="an200_013_ip_vs_lp_saa_benchmarks_300mod5.png", seppngs=separates,\
                        detailpngson=True, polouts=pfiles)
    print "saa benchmarks plotted"

    print "done."