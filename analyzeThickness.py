from glob import glob
from Thickness import Thickness

def analyzeFile(upName, dnName, tileName, frameUp = None, frameDn = None, outdir = "ThicknessOutput"):
    thick = Thickness()

    if (frameUp is not None) and (frameDn is not None):
        thick.loadFrameFiles(frameUp, frameDn)
    thick.loadGelFiles(upName, dnName)
    thick.calculateThickness()
    thick.getSurfaces()
    thick.makeSummary(outdir + "/" + tileName, tileName)
    thick.makeOutput(outdir + "/" + tileName)


if __name__ == "__main__":
    

    upName = "TRIUMF/200220/data/gel_sorted_{tile:d}_up.dat"
    dnName = "TRIUMF/200220/data/gel_sorted_{tile:d}_dn.dat"
    dnName = "TRIUMF/200220/data/gel_sorted_{tile:d}_dn.dat"
    frameUp = "TRIUMF/200220/data/alu-up-{tile:d}.dat"
    frameDn = "TRIUMF/200220/data/alu-dn-{tile:d}.dat"

    for i in range(30):

        try:

            analyzeFile( 
                upName = upName.format(tile = i+1),
                dnName = dnName.format(tile = i+1),
                tileName = f"Tile_{i+1}",
                frameUp = frameUp.format(tile = i+1),
                frameDn = frameDn.format(tile = i+1)
            )
        except FileNotFoundError :
            continue