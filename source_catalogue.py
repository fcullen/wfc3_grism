import figs

def make_source_catalogue():

	if figs.options['FORCE_CATALOG'] is not None:
		### check there's an accompanying segmentation image:
		check_segmentation_image()
		### Assumes that the MAG_AUTO column has already been changed appropriately to something like MAG_F1392W
		sexCat = figs.sex.mySexCat(figs.options['PRE_MADE_INPUT_CATALOGUE'])
		sexCat.write(ROOT_DIRECT+'_drz.cat')
	else:
		pass


def check_segmentation_image():
	if figs.options['PRE_MADE_SEGMENTATION_MAP'] is None:
		figs.showMessage("PRE MADE INPUT CATALOGUE DOES NOT HAVE ACCOMPANYING SEG MAP")
	else:
		pass