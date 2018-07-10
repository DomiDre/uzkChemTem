import os, glob
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import javabridge
import bioformats

def init_vm():
  javabridge.start_vm(class_path=bioformats.JARS)

  #remove annoying logs
  myloglevel="ERROR"  # user string argument for logLevel.
  rootLoggerName = javabridge.get_static_field("org/slf4j/Logger","ROOT_LOGGER_NAME", "Ljava/lang/String;")
  rootLogger = javabridge.static_call("org/slf4j/LoggerFactory","getLogger", "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
  logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level",myloglevel, "Lch/qos/logback/classic/Level;")
  javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)

def kill_vm():
  javabridge.kill_vm()

def load_file(tiffile, pixel_per_nm=None, _init_vm=True, _kill_vm_on_end=True):
  # get metadata
  metadata = os.popen(
    os.path.dirname(os.path.realpath(__file__))+'/showinf ' + tiffile +\
    ' -nopix'
  ).read()

  #Acquisition
  if pixel_per_nm is None:
    if "Nanometers per pixel (X)" in metadata:
      pixel_size_x = float(metadata.split("Nanometers per pixel (X)",1)[-1]\
                    .split('\n')[0]\
                    .split(':')[-1].strip())
      pixel_size_y = float(metadata.split("Nanometers per pixel (Y)",1)[-1]\
                    .split('\n')[0]\
                    .split(':')[-1].strip())
    else:
      print('Unknown TEM file format. Could not read nm per pixel. Set to 1')
      pixel_size_x = 1
      pixel_size_y = 1

    if(pixel_size_x != pixel_size_y):
      print("WARNING. PIXEL SIXE X != PIXEL SIZE Y. UNEXPECTED")
      print(f'Pixel Size X: {pixel_size_x}')
      print(f'Pixel Size Y: {pixel_size_y}')
      print('Continuing with pixel size X only')

    nm_per_pixel = float(pixel_size_x)

  else:
    nm_per_pixel = 1/pixel_per_nm

  # get image
  if _init_vm:
    init_vm()

  #properly select image reader to load data
  #see: https://github.com/CellProfiler/python-bioformats/issues/23
  image_reader = bioformats.formatreader.make_image_reader_class()()
  image_reader.allowOpenToCheckType(True)
  image_reader.setId(tiffile)
  wrapper = bioformats.formatreader.ImageReader(path=tiffile, perform_init=False)
  wrapper.rdr = image_reader
  data = wrapper.read()[::-1,:].T

  if _kill_vm_on_end:
    javabridge.kill_vm()

  Nx, Ny = data.shape

  x = [x*nm_per_pixel for x in range(Nx)]
  y = [x*nm_per_pixel for x in range(Ny)]
  return {
    'x': x,
    'y': y,
    'data': data,
    'nm_per_pixel': nm_per_pixel
  }

def convert_tiffolder(foldername, fileformat='.png', x0=0, y0=0, width=None, height=None, scaleBar=50):
  tiffiles = sorted(glob.glob(foldername + '/*.tif'))
  pretty_plot(tiffiles[0], tiffiles[0].replace('.tif', fileformat), x0, y0, width, height, scaleBar=scaleBar,
    _init_vm=True, _kill_vm_on_end=False, )
  for filename in tiffiles[1:-1]:
    pretty_plot(filename, filename.replace('.tif', fileformat), x0, y0, width, height, scaleBar=scaleBar,
      _init_vm=False, _kill_vm_on_end=False, )
  pretty_plot(tiffiles[-1], tiffiles[-1].replace('.tif', fileformat), x0, y0, width, height, scaleBar=scaleBar,
    _init_vm=False, _kill_vm_on_end=True, )

def pretty_plot(filename, savename=None, x0=0, y0=0, width=None, height=None,
  pixel_per_nm=None, scaleBar=50, _init_vm=True, _kill_vm_on_end=True):
  dataContainer = load_file(filename, pixel_per_nm=pixel_per_nm,
    _init_vm=_init_vm, _kill_vm_on_end=_kill_vm_on_end)
  x = dataContainer['x']
  y = dataContainer['y']
  data = dataContainer['data']
  nm_per_pixel = dataContainer['nm_per_pixel']

  if width is None:
    width = data.shape[0]
  if height is None:
    height = data.shape[1]

  plotX = x[x0:x0+width]
  plotY = y[y0:y0+height]
  plotData = data[x0:x0+width, y0:y0+height]

  fig = plt.figure(frameon=False)
  ax = plt.Axes(fig, [0., 0., 1., 1.], )
  ax.set_axis_off()
  fig.add_axes(ax)
  ax.pcolormesh(plotX, plotY, plotData.T, cmap="gray")
  ax.set_aspect('equal')
  ax.set_xticklabels('')
  ax.set_yticklabels('')
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_xlim(plotX[0], plotX[-1])
  ax.set_ylim(plotY[0], plotY[-1])

  t = ax.text(
    plotX[0]+39.9, plotY[0]+100, str(scaleBar) + ' nm',
    horizontalalignment='center',
    color='black')

  ax.figure.canvas.draw()
  bb = t.get_window_extent(renderer= fig.canvas.renderer)
  transf = ax.transData.inverted()
  transfBb = bb.transformed(transf)
  textWidth = transfBb.width
  textHeight = transfBb.height

  t.remove()
  offsetTextScalebar = max(textWidth, scaleBar)
  rightOffset = 15*nm_per_pixel
  bottomOffset = 15 * nm_per_pixel
  ax.text(
    plotX[-1] - rightOffset - offsetTextScalebar/2,
    plotY[0] + bottomOffset + textHeight*1/4,
    '$'+str(scaleBar)+' \, nm$',\
    horizontalalignment='center',
    color='white')

  ax.add_patch(
    patches.Rectangle(
      (plotX[-1] - rightOffset - offsetTextScalebar/2 - scaleBar*1/2.,\
      plotY[0] + bottomOffset + textHeight*1/20),   # (x,y)
      scaleBar,          # width
      textHeight*1/40,          # height
      color='white'
    )
  )

  if savename is not None:
    save_plot(savename, fig)
  return fig, ax

def save_plot(savename, fig, dpi=None):
  if dpi is not None:
    fig.savefig(savename, dpi=dpi)
  else:
    fig.savefig(savename)
  print("Saved plot to " + savename)