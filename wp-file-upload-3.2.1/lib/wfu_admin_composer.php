<?php

function wfu_shortcode_composer($data = '', $shortcode_tag = 'wordpress_file_upload') {
	global $wpdb;
	global $wp_roles;
	$siteurl = site_url();
 
	$components = wfu_component_definitions();
	if ( $shortcode_tag == 'wordpress_file_upload' ) {
		$cats = wfu_category_definitions();
		$defs = wfu_attribute_definitions();
	}
	else {
		$cats = wfu_browser_category_definitions();
		$defs = wfu_browser_attribute_definitions();
	}
	$plugin_options = wfu_decode_plugin_options(get_option( "wordpress_file_upload_options" ));
	
	if ( $data == "" ) {
		$shortcode = $plugin_options['shortcode'];
		$shortcode_full = '['.$shortcode_tag.' '.$shortcode.']';
		$postid = "";
		$postname = "";
		$posttype = "";
		$posthash = "";
		$shortcode_position = -1;
		$autosave = true;
	}
	else {
		$shortcode = trim(substr($data['shortcode'], strlen('['.$shortcode_tag), -1));
		$shortcode_full = $data['shortcode'];
		$postid = $data['post_id'];
		$postname = get_the_title($postid);
		$posttype = get_post_type($postid);
		$posthash = $data['post_hash'];
		$shortcode_position = $data['position'];
		$autosave = false;
	}
	
	$shortcode_attrs = wfu_shortcode_string_to_array($shortcode);
	$shortcode_id = '';
	foreach ( $defs as $key => $def ) {
		$defs[$key]['default'] = $def['value'];
		if ( array_key_exists($def['attribute'], $shortcode_attrs) ) {
			$defs[$key]['value'] = $shortcode_attrs[$def['attribute']];
		}
		if ( $def['attribute'] == 'uploadid' || $def['attribute'] == 'browserid' ) $shortcode_id = $defs[$key]['value'];
	}

	// index $components
	$components_indexed = array();
	foreach ( $components as $component ) $components_indexed[$component['id']] = $component;
	// index dependiencies
	$governors = array();

	$echo_str = '<div id="wfu_wrapper" class="wrap">';
	$echo_str .= "\n\t".'<h2>Wordpress File Upload Control Panel</h2>';
	$echo_str .= "\n\t".'<div id="wfu_page_obsolete_message" class="error" style="display:none;">';
	$echo_str .= "\n\t\t".'<p>'.WFU_DASHBOARD_PAGE_OBSOLETE.'</p>';
	$echo_str .= "\n\t".'</div>';
	$echo_str .= "\n\t".'<div id="wfu_update_rejected_message" class="error" style="display:none;">';
	$echo_str .= "\n\t\t".'<p>'.WFU_DASHBOARD_UPDATE_SHORTCODE_REJECTED.'</p>';
	$echo_str .= "\n\t".'</div>';
	$echo_str .= "\n\t".'<div id="wfu_update_failed_message" class="error" style="display:none;">';
	$echo_str .= "\n\t\t".'<p>'.WFU_DASHBOARD_UPDATE_SHORTCODE_FAILED.'</p>';
	$echo_str .= "\n\t".'</div>';
	$echo_str .= "\n\t".'<div style="margin-top:20px;">';
	if ( current_user_can( 'manage_options' ) ) $echo_str .= "\n\t".'<a href="'.$siteurl.'/wp-admin/options-general.php?page=wordpress_file_upload&amp;action=manage_mainmenu" class="button" title="go back">Go to Main Menu</a>';
	$echo_str .= "\n\t".'</div>';
	$echo_str .= "\n\t".'<h2 style="margin-bottom: 10px; margin-top: 20px;">Shortcode Composer for '.( $data == "" ? 'Test' : $posttype.' "'.$postname.'" ('.$postid.') Position '.$data['position'] ).' with ID '.$shortcode_id.'</h2>';
	$echo_str .= "\n\t".'<div style="margin-top:10px; display:inline-block;">';
	if ( $data != "") $echo_str .= "\n\t\t".'<input id="wfu_update_shortcode" type="button" value="Update" class="button-primary" disabled="disabled" onclick="wfu_save_shortcode()" /><span id="wfu_update_shortcode_wait" class="spinner" style="float:right; display:none;"></span>';
	$echo_str .= "\n\t\t".'<input id="wfu_shortcode_original_enc" type="hidden" value="'.wfu_plugin_encode_string($shortcode_full).'" />';
	$echo_str .= "\n\t\t".'<input id="wfu_shortcode_tag" type="hidden" value="'.$shortcode_tag.'" />';
	$echo_str .= "\n\t\t".'<input id="wfu_shortcode_postid" type="hidden" value="'.$postid.'" />';
	$echo_str .= "\n\t\t".'<input id="wfu_shortcode_posthash" type="hidden" value="'.$posthash.'" />';
	$echo_str .= "\n\t\t".'<input id="wfu_shortcode_position" type="hidden" value="'.$shortcode_position.'" />';
	$echo_str .= "\n\t".'</div>';
	$echo_str .= "\n\t".'<div style="margin-top:20px;">';
	$echo_str .= "\n\t\t".'<div class="wfu_shortcode_container">';
	$echo_str .= "\n\t\t\t".'<span><strong>Generated Shortcode</strong></span>';
	$echo_str .= "\n\t\t\t".'<span id="wfu_save_label" class="wfu_save_label">saved</span>';
	$echo_str .= "\n\t\t\t".'<textarea id="wfu_shortcode" class="wfu_shortcode" rows="5">['.$shortcode_tag.']</textarea>';
	$echo_str .= "\n\t\t\t".'<div id="wfu_attribute_defaults" style="display:none;">';
	foreach ( $defs as $def )
		$echo_str .= "\n\t\t\t\t".'<input id="wfu_attribute_default_'.$def['attribute'].'" type="hidden" value="'.$def['default'].'" />';
	$echo_str .= "\n\t\t\t".'</div>';
	$echo_str .= "\n\t\t\t".'<div id="wfu_attribute_values" style="display:none;">';
	foreach ( $defs as $def )
		$echo_str .= "\n\t\t\t\t".'<input id="wfu_attribute_value_'.$def['attribute'].'" type="hidden" value="'.$def['value'].'" />';
	$echo_str .= "\n\t\t\t".'</div>';
	$echo_str .= "\n\t\t".'</div>';
	$echo_str .= "\n\t".'</div>';
	$echo_str .= "\n\t".'<h3 id="wfu_tab_container" class="nav-tab-wrapper">';
	$is_first = true;
	foreach ( $cats as $key => $cat ) {
		$echo_str .= "\n\t\t".'<a id="wfu_tab_'.$key.'" class="nav-tab'.( $is_first ? ' nav-tab-active' : '' ).'" href="javascript: wfu_admin_activate_tab(\''.$key.'\');">'.$cat.'</a>';
		$is_first = false;
	}
	$echo_str .= "\n\t".'</h3>';

	$prevcat = "";
	$prevsubcat = "";
	$is_first = true;
	$block_open = false;
	$subblock_open = false;
	foreach ( $defs as $def ) {
		$attr = $def['attribute'];
		$subblock_active = false;
		//detect if the dependencies of this attribute will be disabled or not
		if ( ( $def['type'] == "onoff" && $def['value'] == "true" ) ||
			( $def['type'] == "radio" && in_array("*".$def['value'], $def['listitems']) ) )
			$subblock_active = true;
		// assign dependencies if exist
		if ( $def['dependencies'] != null )
			foreach ( $def['dependencies'] as $dependency ) {
				if ( substr($dependency, 0, 1) == "!" ) //invert state for this dependency if an exclamation mark is defined
					$governors[substr($dependency, 1)] = array( 'attribute' => $attr, 'active' => !$subblock_active, 'inv' => '_inv' );
				else
					$governors[$dependency] = array( 'attribute' => $attr, 'active' => $subblock_active, 'inv' => '' );
			}
		//check if this attribute depends on other
		if ( !array_key_exists($attr, $governors) ) $governors[$attr] = "";
		if ( $governors[$attr] != "" ) $governor = $governors[$attr];
		else $governor = array( 'attribute' => "independent", 'active' => true, 'inv' => '' );

		//close previous blocks
		if ( $def['parent'] == "" ) {
			if ( $subblock_open ) {
				$echo_str .= "\n\t\t\t\t\t\t\t".'</tbody>';
				$echo_str .= "\n\t\t\t\t\t\t".'</table>';
				$subblock_open = false;
			}
			if ( $block_open ) {
				$echo_str .= "\n\t\t\t\t\t".'</div></td>';
				$echo_str .= "\n\t\t\t\t".'</tr>';
				$block_open = false;
			}
		}
		//check if new category must be generated
		if ( $def['category'] != $prevcat ) {
			if ( $prevcat != "" ) {
				$echo_str .= "\n\t\t\t".'</tbody>';
				$echo_str .= "\n\t\t".'</table>';
				$echo_str .= "\n\t".'</div>';
			}
			$prevcat = $def['category'];
			$prevsubcat = "";
			$echo_str .= "\n\t".'<div id="wfu_container_'.$prevcat.'" class="wfu_container"'.( $is_first ? '' : ' style="display:none;"' ).'">';
			$echo_str .= "\n\t\t".'<table class="form-table wfu_main_table">';
			$echo_str .= "\n\t\t\t".'<thead><tr><th></th><td></td><td></td></tr></thead>';
			$echo_str .= "\n\t\t\t".'<tbody>';
			$is_first = false;
		}
		//check if new sub-category must be generated
		if ( $def['subcategory'] != $prevsubcat ) {
			$prevsubcat = $def['subcategory'];
			$echo_str .= "\n\t\t\t\t".'<tr class="wfu_subcategory">';
			$echo_str .= "\n\t\t\t\t\t".'<th scope="row" colspan="3">';
			$echo_str .= "\n\t\t\t\t\t\t".'<h3 style="margin-bottom: 10px; margin-top: 10px;">'.$prevsubcat.'</h3>';
			$echo_str .= "\n\t\t\t\t\t".'</th>';
			$echo_str .= "\n\t\t\t\t".'</tr>';
		}
		//draw attribute element
		if ( $def['parent'] == "" ) {
			$dlp = "\n\t\t\t\t";
		}
		else {
			if ( !$subblock_open ) {
				$echo_str .= "\n\t\t\t\t\t\t".'<div class="wfu_shadow wfu_shadow_'.$def['parent'].$governor['inv'].'" style="display:'.( $governor['active'] ? 'none' : 'block' ).';"></div>';
				$echo_str .= "\n\t\t\t\t\t\t".'<table class="form-table wfu_inner_table" style="margin:0;">';
				$echo_str .= "\n\t\t\t\t\t\t\t".'<tbody>';
			}
			$dlp = "\n\t\t\t\t\t\t\t\t";
		}
		$echo_str .= $dlp.'<tr>';
		$echo_str .= $dlp."\t".'<th scope="row"><div class="wfu_td_div">';
		if ( $def['parent'] == "" ) $echo_str .= $dlp."\t\t".'<div class="wfu_shadow wfu_shadow_'.$governor['attribute'].$governor['inv'].'" style="display:'.( $governor['active'] ? 'none' : 'block' ).';"></div>';
		$echo_str .= $dlp."\t\t".'<div class="wfu_restore_container" title="Double-click to restore defaults setting"><img src="'.WFU_IMAGE_ADMIN_RESTOREDEFAULT.'" ondblclick="wfu_apply_value(\''.$attr.'\', \''.$def['type'].'\', \''.$def['default'].'\');" /></div>';
		$echo_str .= $dlp."\t\t".'<label for="wfu_attribute_'.$attr.'">'.$def['name'].'</label>';
		$echo_str .= $dlp."\t\t".'<div class="wfu_help_container" title="'.$def['help'].'"><img src="'.WFU_IMAGE_ADMIN_HELP.'" /></div>';
		$echo_str .= $dlp."\t".'</div></th>';
		$echo_str .= $dlp."\t".'<td style="vertical-align:top;"><div class="wfu_td_div">';
		if ( $def['parent'] == "" ) $echo_str .= $dlp."\t\t".'<div class="wfu_shadow wfu_shadow_'.$governor['attribute'].$governor['inv'].'" style="display:'.( $governor['active'] ? 'none' : 'block' ).';"></div>';
		if ( $def['type'] == "onoff" ) {
			$echo_str .= $dlp."\t\t".'<div id="wfu_attribute_'.$attr.'" class="wfu_onoff_container_'.( $def['value'] == "true" ? "on" : "off" ).'" onclick="wfu_admin_onoff_clicked(\''.$attr.'\');">';
			$echo_str .= $dlp."\t\t\t".'<div class="wfu_onoff_slider"></div>';
			$echo_str .= $dlp."\t\t\t".'<span class="wfu_onoff_text">ON</span>';
			$echo_str .= $dlp."\t\t\t".'<span class="wfu_onoff_text">OFF</span>';
			$echo_str .= $dlp."\t\t".'</div>';
		}
		elseif ( $def['type'] == "text" ) {
			$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
			$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="text" name="wfu_text_elements" value="'.$val.'" />';
			if ( $def['variables'] != null ) $echo_str .= $dlp.wfu_insert_variables($def['variables'], 'wfu_variable wfu_variable_'.$attr);
		}
		elseif ( $def['type'] == "placements" ) {
			$components_used = array();
			foreach ( $components as $component ) $components_used[$component['id']] = false;
			$centered_content = '<div style="display:table; width:100%; height:100%;"><div style="display:table-cell; text-align:center; vertical-align:middle;">XXX</div></div>';
			$echo_str .= $dlp."\t\t".'<div class="wfu_placements_wrapper">';
			$echo_str .= $dlp."\t\t\t".'<div id="wfu_placements_container" class="wfu_placements_container">';
			$itemplaces = explode("/", $def['value']);
			foreach ( $itemplaces as $section ) {
				$echo_str .= $dlp."\t\t\t\t".'<div class="wfu_component_separator_hor"></div>';
				$echo_str .= $dlp."\t\t\t\t".'<div class="wfu_component_separator_ver"></div>';
				$items_in_section = explode("+", trim($section));
				$section_array = array( );
				foreach ( $items_in_section as $item_in_section ) {
					if ( key_exists($item_in_section, $components_indexed) ) {
						$components_used[$item_in_section] = true;
						$echo_str .= $dlp."\t\t\t\t".'<div id="wfu_component_box_'.$item_in_section.'" class="wfu_component_box" draggable="true" title="'.$components_indexed[$item_in_section]['help'].'">'.str_replace("XXX", $components_indexed[$item_in_section]['name'], $centered_content).'</div>';
						$echo_str .= $dlp."\t\t\t\t".'<div class="wfu_component_separator_ver"></div>';
					}
				}
			}
			$echo_str .= $dlp."\t\t\t\t".'<div class="wfu_component_separator_hor"></div>';
			$echo_str .= $dlp."\t\t\t\t".'<div id="wfu_component_bar_hor" class="wfu_component_bar_hor"></div>';
			$echo_str .= $dlp."\t\t\t\t".'<div id="wfu_component_bar_ver" class="wfu_component_bar_ver"></div>';
			$echo_str .= $dlp."\t\t\t".'</div>';
			$echo_str .= $dlp."\t\t\t".'<div id="wfu_componentlist_container" class="wfu_componentlist_container">';
			$echo_str .= $dlp."\t\t\t\t".'<div id="wfu_componentlist_dragdrop" class="wfu_componentlist_dragdrop" style="display:none;"></div>';
			$ii = 1;
			foreach ( $components as $component ) {
				$echo_str .= $dlp."\t\t\t\t".'<div id="wfu_component_box_container_'.$component['id'].'" class="wfu_component_box_container">';
				$echo_str .= $dlp."\t\t\t\t\t".'<div class="wfu_component_box_base">'.str_replace("XXX", $component['name'], $centered_content).'</div>';
				if ( !$components_used[$component['id']] )
					$echo_str .= $dlp."\t\t\t\t\t".'<div id="wfu_component_box_'.$component['id'].'" class="wfu_component_box wfu_inbase" draggable="true" title="'.$component['help'].'">'.str_replace("XXX", $component['name'], $centered_content).'</div>';
				$echo_str .= $dlp."\t\t\t\t".'</div>'.( ($ii++) % 3 == 0 ? '<br />' : '' );
			}
			$echo_str .= $dlp."\t\t\t".'</div>';
			$echo_str .= $dlp."\t\t".'</div>';
		}
		elseif ( $def['type'] == "ltext" ) {
			$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
			$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="text" name="wfu_text_elements" class="wfu_long_text" value="'.$val.'" />';
			if ( $def['variables'] != null ) $echo_str .= $dlp.wfu_insert_variables($def['variables'], 'wfu_variable wfu_variable_'.$attr);
		}
		elseif ( $def['type'] == "integer" ) {
			$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
			$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="number" name="wfu_text_elements" class="wfu_short_text" min="1" value="'.$val.'" />';
			if ( isset($def['listitems']['unit']) ) $echo_str .= $dlp."\t\t".'<label> '.$def['listitems']['unit'].'</label>';
		}
		elseif ( $def['type'] == "float" ) {
			$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
			$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="number" name="wfu_text_elements" class="wfu_short_text" step="any" min="0" value="'.$val.'" />';
			if ( isset($def['listitems']['unit']) ) $echo_str .= $dlp."\t\t".'<label> '.$def['listitems']['unit'].'</label>';
		}
		elseif ( $def['type'] == "date" ) {
			$val = $def['value'];
			$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="date" name="wfu_date_elements" value="'.$val.'" />';
		}
		elseif ( $def['type'] == "radio" ) {
			$echo_str .= $dlp."\t\t";
			$ii = 0;
			foreach ( $def['listitems'] as $item )
				$echo_str .= '<input name="wfu_radioattribute_'.$attr.'" type="radio" value="'.$item.'" '.( $item == $def['value'] || $item == "*".$def['value'] ? 'checked="checked" ' : '' ).'style="width:auto; margin:0px 2px 0px '.( ($ii++) == 0 ? '0px' : '8px' ).';" onchange="wfu_admin_radio_clicked(\''.$attr.'\');" />'.( $item[0] == "*" ? substr($item, 1) : $item );
//			$echo_str .= '<input type="button" class="button" value="empty" style="width:auto; margin:-2px 0px 0px 8px;" />';
		}
		elseif ( $def['type'] == "ptext" ) {
			$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
			$parts = explode("/", $val);
			$singular = $parts[0];
			if ( count($parts) < 2 ) $plural = $singular;
			else $plural = $parts[1];
			$echo_str .= $dlp."\t\t".'<span class="wfu_ptext_span">Singular</span><input id="wfu_attribute_s_'.$attr.'" type="text" name="wfu_ptext_elements" value="'.$singular.'" />';
			if ( $def['variables'] != null ) if ( count($def['variables']) > 0 ) $echo_str .= $dlp."\t\t".'<br /><span class="wfu_ptext_span">&nbsp;</span>';
			if ( $def['variables'] != null ) $echo_str .= $dlp.wfu_insert_variables($def['variables'], 'wfu_variable wfu_variable_s_'.$attr);
			$echo_str .= $dlp."\t\t".'<br /><span class="wfu_ptext_span">Plural</span><input id="wfu_attribute_p_'.$attr.'" type="text" name="wfu_ptext_elements" value="'.$plural.'" />';
			if ( $def['variables'] != null ) if ( count($def['variables']) > 0 ) $echo_str .= $dlp."\t\t".'<br /><span class="wfu_ptext_span">&nbsp;</span>';
			if ( $def['variables'] != null ) $echo_str .= $dlp.wfu_insert_variables($def['variables'], 'wfu_variable wfu_variable_p_'.$attr, $dlp);
		}
		elseif ( $def['type'] == "mtext" ) {
			$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
			$echo_str .= $dlp."\t\t".'<textarea id="wfu_attribute_'.$attr.'" name="wfu_text_elements" rows="5">'.$val.'</textarea>';
			if ( $def['variables'] != null ) $echo_str .= $dlp.wfu_insert_variables($def['variables'], 'wfu_variable wfu_variable_'.$attr);
		}
		elseif ( $def['type'] == "folderlist" ) {
			$echo_str .= $dlp."\t\t".'<div id="wfu_subfolders_inner_shadow_'.$attr.'" class="wfu_subfolders_inner_shadow" style="display:none;"></div>';
			$subfolders = wfu_parse_folderlist($def['value']);
			$poptitle = "Populate list automatically with the first-level subfolders of the path defined in uploadpath";
			$edittitle = "Allow the user to type the subfolder and filter the list during typing";
			$echo_str .= $dlp."\t\t".'<input type="checkbox" id="wfu_subfolders_auto_'.$attr.'"'.( substr($def['value'], 0, 4) == "auto" ? ' checked="checked"' : '' ).' onchange="wfu_subfolders_auto_changed(\''.$attr.'\');" title="'.$poptitle.'" /><label for="wfu_subfolders_auto_'.$attr.'" title="'.$poptitle.'"> Auto-populate list</label>';
			$echo_str .= $dlp."\t\t".'<div style="display:'.( substr($def['value'], 0, 4) == "auto" ? 'inline' : 'none' ).'; padding:0; margin:0 0 0 30px; background:none; border:none;"><input type="checkbox" id="wfu_subfolders_editable_'.$attr.'"'.( substr($def['value'], 0, 5) == "auto+" ? ' checked="checked"' : '' ).' onchange="wfu_subfolders_auto_changed(\''.$attr.'\');" title="'.$edittitle.'" /><label for="wfu_subfolders_editable_'.$attr.'" title="'.$edittitle.'"> List is editable</label></div><br />';
			$echo_str .= $dlp."\t\t".'<input type="hidden" id="wfu_subfolders_manualtext_'.$attr.'" value="'.( substr($def['value'], 0, 4) == "auto" ? "" : $def['value'] ).'" />';
			$echo_str .= $dlp."\t\t".'<select id="wfu_attribute_'.$attr.'" class="wfu_select_folders'.( count($subfolders['path']) == 0 ? ' wfu_select_folders_empty' : '' ).'" size="7"'.( substr($def['value'], 0, 4) == "auto" ? ' disabled="disabled"' : '' ).' onchange="wfu_subfolders_changed(\''.$attr.'\');">';
			foreach ($subfolders['path'] as $ind => $subfolder) {
				if ( substr($subfolder, -1) == '/' ) $subfolder = substr($subfolder, 0, -1);
				$subfolder_raw = explode('/', $subfolder);
				$subfolder = $subfolder_raw[count($subfolder_raw) - 1];
				$text = str_repeat("&nbsp;&nbsp;&nbsp;", intval($subfolders['level'][$ind])).$subfolders['label'][$ind];
				$subvalue = str_repeat("*", intval($subfolders['level'][$ind])).( $subfolders['default'][$ind] ? '&' : '' ).( $subfolder == "" ? '{root}' : $subfolder ).'/'.$subfolders['label'][$ind];
				$echo_str .= $dlp."\t\t\t".'<option class="'.( $subfolders['default'][$ind] ? 'wfu_select_folders_option_default' : '' ).'" value="'.wfu_plugin_encode_string($subvalue).'">'.$text.'</option>';
			}
			$echo_str .= $dlp."\t\t\t".'<option value="">'.( substr($def['value'], 0, 4) != "auto" && count($subfolders['path']) == 0 ? 'press here' : '' ).'</option>';
			$echo_str .= $dlp."\t\t".'</select>';
			$echo_str .= $dlp."\t\t".'<div id="wfu_subfolder_nav_'.$attr.'" class="wfu_subfolder_nav_container">';
			$echo_str .= $dlp."\t\t\t".'<table class="wfu_subfolder_nav"><tbody>';
			$echo_str .= $dlp."\t\t\t\t".'<tr><td><button id="wfu_subfolders_up_'.$attr.'" name="wfu_subfolder_nav_'.$attr.'" class="button" disabled="disabled" title="move item up" onclick="wfu_subfolders_up_clicked(\''.$attr.'\');">&uarr;</button></tr></td>';
			$echo_str .= $dlp."\t\t\t\t".'<tr><td><button id="wfu_subfolders_left_'.$attr.'" name="wfu_subfolder_nav_'.$attr.'" class="button" title="make it parent" disabled="disabled" style="height:14px;" onclick="wfu_subfolders_left_clicked(\''.$attr.'\');">&larr;</button>';
			$echo_str .= $dlp."\t\t\t\t".'<button id="wfu_subfolders_right_'.$attr.'" name="wfu_subfolder_nav_'.$attr.'" class="button" title="make it child" disabled="disabled" style="height:14px;" onclick="wfu_subfolders_right_clicked(\''.$attr.'\');">&rarr;</button></tr></td>';
			$echo_str .= $dlp."\t\t\t\t".'<tr><td><button id="wfu_subfolders_down_'.$attr.'" name="wfu_subfolder_nav_'.$attr.'" class="button" title="move item down" disabled="disabled" onclick="wfu_subfolders_down_clicked(\''.$attr.'\');">&darr;</button></tr></td>';
			$echo_str .= $dlp."\t\t\t\t".'<tr><td style="line-height:0;"><button  class="button" style="visibility:hidden; height:10px;"></button></tr></td>';
			$echo_str .= $dlp."\t\t\t\t".'<tr><td><button id="wfu_subfolders_add_'.$attr.'" name="wfu_subfolder_nav_'.$attr.'" class="button" title="add new item" disabled="disabled" style="height:14px;" onclick="wfu_subfolders_add_clicked(\''.$attr.'\');">+</button></tr></td>';
			$echo_str .= $dlp."\t\t\t\t".'<tr><td><button id="wfu_subfolders_def_'.$attr.'" name="wfu_subfolder_nav_'.$attr.'" class="button" title="make it default" disabled="disabled" style="height:14px;" onclick="wfu_subfolders_def_clicked(\''.$attr.'\');">&diams;</button></tr></td>';
			$echo_str .= $dlp."\t\t\t\t".'<tr><td><button id="wfu_subfolders_del_'.$attr.'" name="wfu_subfolder_nav_'.$attr.'" class="button" title="delete item" disabled="disabled" style="height:14px;" onclick="wfu_subfolders_del_clicked(\''.$attr.'\');">-</button></tr></td>';
			$echo_str .= $dlp."\t\t\t".'</tbody></table>';
			$echo_str .= $dlp."\t\t".'</div>';
			$echo_str .= $dlp."\t\t".'<div id="wfu_subfolder_tools_'.$attr.'" class="wfu_subfolder_tools_container wfu_subfolder_tools_disabled">';
			$echo_str .= $dlp."\t\t\t".'<table class="wfu_subfolder_tools"><tbody><tr>';
			$echo_str .= $dlp."\t\t\t\t".'<td style="width:40%;">';
			$echo_str .= $dlp."\t\t\t\t\t".'<label>Label</label>';
			$echo_str .= $dlp."\t\t\t\t\t".'<input id="wfu_subfolders_label_'.$attr.'" name="wfu_subfolder_tools_input" type="text" disabled="disabled" />';
			$echo_str .= $dlp."\t\t\t\t".'</td>';
			$echo_str .= $dlp."\t\t\t\t".'<td style="width:60%;"><div style="padding-right:36px;">';
			$echo_str .= $dlp."\t\t\t\t\t".'<label>Path</label>';
			$echo_str .= $dlp."\t\t\t\t\t".'<input id="wfu_subfolders_path_'.$attr.'" name="wfu_subfolder_tools_input" type="text" disabled="disabled" />';
			$echo_str .= $dlp."\t\t\t\t\t".'<button id="wfu_subfolders_browse_'.$attr.'" class="button" title="browse folders" style="right:18px;" disabled="disabled" onclick="wfu_subfolders_browse_clicked(\''.$attr.'\');"><img src="'.WFU_IMAGE_ADMIN_SUBFOLDER_BROWSE.'" ></button>';
			$echo_str .= $dlp."\t\t\t\t\t".'<button id="wfu_subfolders_ok_'.$attr.'" class="button" title="save changes" style="right:0px;" disabled="disabled" onclick="wfu_subfolders_ok_clicked(\''.$attr.'\');"><img src="'.WFU_IMAGE_ADMIN_SUBFOLDER_OK.'" ></button>';
			// file browser dialog
			$echo_str .= $dlp."\t\t\t\t\t".'<div id="wfu_subfolders_browser_'.$attr.'" class="wfu_subfolders_browser_container" style="display:none;">';
			$echo_str .= $dlp."\t\t\t\t\t\t".'<table><tbody>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t".'<tr><td style="height:15px;">';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t".'<div>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t".'<label>Folder Browser</label>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t".'<button class="button wfu_folder_browser_cancel" onclick="wfu_folder_browser_cancel_clicked(\''.$attr.'\');"><img src="'.WFU_IMAGE_ADMIN_SUBFOLDER_CANCEL.'" ></button>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t".'</div>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t".'</td></tr>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t".'<tr><td style="height:106px;">';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t".'<div>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t".'<select id="wfu_subfolders_browser_list_'.$attr.'" size="2" onchange="wfu_subfolders_browser_list_changed(\''.$attr.'\');">';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t\t".'<option>Value</option>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t\t".'<option>Value2</option>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t\t".'<option>Value3</option>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t".'</select>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t".'<div id="wfu_subfolders_browser_msgcont_'.$attr.'" class="wfu_folder_browser_loading_container" style="padding-top:40px;">';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t\t".'<label id="wfu_subfolders_browser_msg_'.$attr.'" style="margin-bottom:4px;">loading folder contents...</label>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t\t".'<img id="wfu_subfolders_browser_img_'.$attr.'" src="'.WFU_IMAGE_ADMIN_SUBFOLDER_LOADING.'" ></button>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t".'</div>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t".'</div>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t".'</td></tr>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t".'<tr><td align="right" style="height:15px;">';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t".'<div>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t".'<button class="button" onclick="wfu_folder_browser_cancel_clicked(\''.$attr.'\');">Cancel</button>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t\t".'<button id="wfu_subfolders_browser_ok_'.$attr.'" class="button">Ok</button>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t\t".'</div>';
			$echo_str .= $dlp."\t\t\t\t\t\t\t".'</td></tr>';
			$echo_str .= $dlp."\t\t\t\t\t\t".'</tbody></table>';
			$echo_str .= $dlp."\t\t\t\t\t".'</div>';

			$echo_str .= $dlp."\t\t\t\t".'</div></td>';
			$echo_str .= $dlp."\t\t\t".'</tr></tbody></table>';
			$echo_str .= $dlp."\t\t\t".'<input id="wfu_subfolders_isnewitem_'.$attr.'" type="hidden" value="" />';
			$echo_str .= $dlp."\t\t\t".'<input id="wfu_subfolders_newitemindex_'.$attr.'" type="hidden" value="" />';
			$echo_str .= $dlp."\t\t\t".'<input id="wfu_subfolders_newitemlevel_'.$attr.'" type="hidden" value="" />';
			$echo_str .= $dlp."\t\t\t".'<input id="wfu_subfolders_newitemlevel2_'.$attr.'" type="hidden" value="" />';
			$echo_str .= $dlp."\t\t".'</div>';
		}
		elseif ( $def['type'] == "mchecklist" ) {
			$help_count = 0;
			foreach ( $def['listitems'] as $key => $item ) {
				$parts = explode("/", $item);
				if ( count($parts) == 1 ) {
					$items[$key]['id'] = $item;
					$items[$key]['help'] = '';
				}
				else {
					$items[$key]['id'] = $parts[0];
					$items[$key]['help'] = $parts[1];
					$help_count ++;
				}
			}
			$def['value'] = strtolower($def['value']);
			if ( $def['value'] == "all" ) $selected = array();
			else $selected = explode(",", $def['value']);
			foreach ( $selected as $key => $item ) $selected[$key] = trim($item);
			$echo_str .= $dlp."\t\t".'<div id="wfu_attribute_'.$attr.'" class="wfu_mchecklist_container">';
			$is_first = true;
			foreach ( $items as $key => $item ) {
				if ( !$is_first ) $echo_str .= "<br />";
				$is_first = false;
				$echo_str .= $dlp."\t\t\t".'<div class="wfu_mchecklist_item"><input id="wfu_attribute_'.$attr.'_'.$key.'" type="checkbox"'.( $def['value'] == "all" || in_array($item['id'], $selected) ? ' checked="checked"' : '' ).( $def['value'] == "all" ? ' disabled="disabled"' : '' ).' onchange="wfu_update_mchecklist_value(\''.$attr.'\');" /><label for="wfu_attribute_'.$attr.'_'.$key.'">'.$item['id'].'</label>';
				if ( $item['help'] != '' ) $echo_str .= '<div class="wfu_help_container" title="'.$item['help'].'"><img src="'.WFU_IMAGE_ADMIN_HELP.'" /></div>';
				$echo_str .= '</div>';
			}
			$echo_str .= $dlp."\t\t".'</div>';
			$echo_str .= $dlp."\t\t".'<div id="wfu_attribute_'.$attr.'_optionhelp" class="wfu_help_container" title="" style="display:none; position:absolute;"><img src="'.WFU_IMAGE_ADMIN_HELP.'" style="visibility:visible;" /></div>';
			$echo_str .= $dlp."\t\t".'<div class="wfu_mchecklist_checkall"><input id="wfu_attribute_'.$attr.'_all" type="checkbox" onchange="wfu_update_mchecklist_value(\''.$attr.'\');"'.( $def['value'] == "all" ? ' checked="checked"' : '' ).' /> Select all</div>';
		}
		elseif ( $def['type'] == "rolelist" ) {
			$roles = $wp_roles->get_names();
			$selected = explode(",", $def['value']);
			$default_administrator = ( is_array($def['listitems']) && in_array('default_administrator', $def['listitems']) );
			if ( in_array('all', $selected) ) $rolesselected = ( $default_administrator ? array("administrator") : array( ) );
			else $rolesselected = $selected;
			foreach ( $selected as $key => $item ) $selected[$key] = trim($item);
			$echo_str .= $dlp."\t\t".'<table class="wfu_rolelist_container"><tbody><tr><td>';
			$echo_str .= $dlp."\t\t".'<select id="wfu_attribute_'.$attr.'" multiple="multiple" size="'.count($roles).'" onchange="wfu_update_rolelist_value(\''.$attr.'\');"'.( in_array('all', $selected) ? ' disabled="disabled"' : '' ).'>';
			foreach ( $roles as $roleid => $rolename )
				$echo_str .= $dlp."\t\t\t".'<option value="'.$roleid.'"'.( in_array($roleid, $rolesselected) ? ' selected="selected"' : '' ).'>'.$rolename.'</option>';
			$echo_str .= $dlp."\t\t".'</select>';
			$echo_str .= $dlp."\t\t".'</td><td>';
			$echo_str .= $dlp."\t\t".'<div class="wfu_rolelist_checkbtn"><input class="'.( $default_administrator ? 'wfu_default_administrator' : '' ).'" id="wfu_attribute_'.$attr.'_all" type="checkbox" onchange="wfu_update_rolelist_value(\''.$attr.'\');"'.( in_array('all', $selected) ? ' checked="checked"' : '' ).' /><label for="wfu_attribute_'.$attr.'_all"> Select all</label></div><br />';
			$echo_str .= $dlp."\t\t".'<div class="wfu_rolelist_checkbtn"><input id="wfu_attribute_'.$attr.'_guests" type="checkbox" onchange="wfu_update_rolelist_value(\''.$attr.'\');"'.( in_array("guests", $selected) ? ' checked="checked"' : '' ).' /><label for="wfu_attribute_'.$attr.'_guests"> Include guests</label></div>';
			$echo_str .= $dlp."\t\t".'</td></tr></tbody></table>';
		}
		elseif ( $def['type'] == "userlist" ) {
			$users = get_users();
			$selected = explode(",", $def['value']);
			$default_0 = ( is_array($def['listitems']) && in_array('default_0', $def['listitems']) );
			if ( in_array('all', $selected) ) $usersselected = ( $default_0 ? array($users[0]->user_login) : array( ) );
			else $usersselected = $selected;
			$only_current = false;
			$echo_str .= $dlp."\t\t".'<table class="wfu_userlist_container"><tbody><tr>';
			if ( is_array($def['listitems']) && in_array('include_current', $def['listitems']) ) {
				$only_current = ( $def['value'] == 'current' );
				if ( $only_current ) $usersselected = ( $default_0 ? array($users[0]->user_login) : array( ) );
				$echo_str .= $dlp."\t\t".'<td colspan="2"><div class="wfu_userlist_checkbtn"><input id="wfu_attribute_'.$attr.'_current" type="checkbox" onchange="wfu_update_userlist_value(\''.$attr.'\');"'.( $only_current ? ' checked="checked"' : '' ).' /><label for="wfu_attribute_'.$attr.'_current"> Only From Current User</label></div>';
				$echo_str .= $dlp."\t\t".'</td></tr><tr>';
			}
			$echo_str .= $dlp."\t\t".'<td><select id="wfu_attribute_'.$attr.'" multiple="multiple" size="'.min(count($users), 10).'" onchange="wfu_update_userlist_value(\''.$attr.'\');"'.( $only_current || in_array('all', $selected) ? ' disabled="disabled"' : '' ).'>';
			foreach ( $users as $userid => $user )
				$echo_str .= $dlp."\t\t\t".'<option value="'.$user->user_login.'"'.( in_array($user->user_login, $usersselected) ? ' selected="selected"' : '' ).'>'.$user->display_name.' ('.$user->user_login.')</option>';
			$echo_str .= $dlp."\t\t".'</select>';
			$echo_str .= $dlp."\t\t".'</td><td>';
			$echo_str .= $dlp."\t\t".'<div class="wfu_userlist_checkbtn"><input class="'.( $default_0 ? 'wfu_default_0' : '' ).'" id="wfu_attribute_'.$attr.'_all" type="checkbox" onchange="wfu_update_userlist_value(\''.$attr.'\');"'.( in_array('all', $selected) ? ' checked="checked"' : '' ).( $only_current ? ' disabled="disabled"' : '' ).' /><label for="wfu_attribute_'.$attr.'_all"> Select all</label></div><br />';
			$echo_str .= $dlp."\t\t".'<div class="wfu_userlist_checkbtn"><input id="wfu_attribute_'.$attr.'_guests" type="checkbox" onchange="wfu_update_userlist_value(\''.$attr.'\');"'.( in_array("guests", $selected) ? ' checked="checked"' : '' ).( $only_current ? ' disabled="disabled"' : '' ).' /><label for="wfu_attribute_'.$attr.'_guests"> Include guests</label></div>';
			$echo_str .= $dlp."\t\t".'</td></tr></tbody></table>';
		}
		elseif ( $def['type'] == "postlist" ) {
			$processed = false;
			if ( is_array($def['listitems']) ) {
				$has_current = in_array('include_current', $def['listitems']);
				if ( $has_current ) unset($def['listitems'][array_search('include_current', $def['listitems'])]);
				foreach ( $def['listitems'] as $post_type ) {
					// if a post type cannot be found then we reset the list so that it is not processed at all
					if ( get_post_type_object( $post_type ) == null ) {
						$def['listitems'] = array();
						break;
					}
				}
				if ( count($def['listitems']) > 0 ) {
					$selected = explode(",", $def['value']);
					$only_current = false;
					$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'_postlist" type="hidden" value="'.implode(",", $def['listitems']).'" />';
					$echo_str .= $dlp."\t\t".'<table class="wfu_postlist_container"><tbody><tr>';
					if ( $has_current ) {
						$only_current = ( $def['value'] == 'current' );
						if ( $only_current ) $sselected = array();
						$echo_str .= $dlp."\t\t".'<td colspan="'.count($def['listitems']).'"><div class="wfu_postlist_checkbtn"><input id="wfu_attribute_'.$attr.'_current" type="checkbox" onchange="wfu_update_postlist_value(\''.$attr.'\');"'.( $only_current ? ' checked="checked"' : '' ).' /><label for="wfu_attribute_'.$attr.'_current"> Only From Current Post/Page</label></div>';
						$echo_str .= $dlp."\t\t".'</td></tr><tr>';
					}
					$postargs = array( 'post_type' => $def['listitems'], 'post_status' => "publish,private,draft", 'posts_per_page' => -1 );
					$posts = get_posts($postargs);
					$list = wfu_construct_post_list($posts);
					foreach ( $def['listitems'] as $post_type ) {
						$flatlist = wfu_flatten_post_list($list[$post_type]);
						$postobj = get_post_type_object( $post_type );
						$echo_str .= $dlp."\t\t".'<td><div class="wfu_postlist_header"><label>'.$postobj->label.'</label><div class="wfu_postlist_selectall"><input id="wfu_attribute_'.$attr.'_all_'.$post_type.'" type="checkbox" onchange="wfu_update_postlist_value(\''.$attr.'\');"'.( in_array('all', $selected) || in_array('all'.$post_type, $selected) ? ' checked="checked"' : '' ).( $only_current ? ' disabled="disabled"' : '' ).' /><label for="wfu_attribute_'.$attr.'_all_'.$post_type.'"> Select all</label></div></div>';
						$echo_str .= $dlp."\t\t".'<select id="wfu_attribute_'.$attr.'_'.$post_type.'" multiple="multiple" size="'.min(count($flatlist), 10).'" onchange="wfu_update_postlist_value(\''.$attr.'\');"'.( $only_current || in_array('all', $selected) || in_array('all'.$post_type, $selected) ? ' disabled="disabled"' : '' ).'>';
						foreach ( $flatlist as $item )
							$echo_str .= $dlp."\t\t\t".'<option value="'.$item['id'].'"'.( in_array($item['id'], $selected) ? ' selected="selected"' : '' ).'>'.str_repeat('&nbsp;', 4 * $item['level']).( $item['status'] == 1 ? '[Private]' : ( $item['status'] == 2 ? '[Draft]' : '' ) ).$item['title'].'</option>';
						$echo_str .= $dlp."\t\t".'</select></td>';
					}
					$echo_str .= $dlp."\t\t".'</tr></tbody></table>';
					$processed = true;
				}
			}
			if ( !processed ) {
				$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
				$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="text" name="wfu_text_elements" value="'.$val.'" />';				
			}
		}
		elseif ( $def['type'] == "bloglist" ) {
			if ( function_exists('wp_get_sites') ) {
				$blogs = wp_get_sites( );
				$selected = explode(",", $def['value']);
				if ( in_array('all', $selected) ) $blogsselected = array( );
				else $blogsselected = $selected;
				$only_current = false;
				$echo_str .= $dlp."\t\t".'<table class="wfu_bloglist_container"><tbody><tr>';
				if ( is_array($def['listitems']) && in_array('include_current', $def['listitems']) ) {
					$only_current = ( $def['value'] == 'current' );
					if ( $only_current ) $blogsselected = array( );
					$echo_str .= $dlp."\t\t".'<td colspan="2"><div class="wfu_bloglist_checkbtn"><input id="wfu_attribute_'.$attr.'_current" type="checkbox" onchange="wfu_update_bloglist_value(\''.$attr.'\');"'.( $only_current ? ' checked="checked"' : '' ).' /><label for="wfu_attribute_'.$attr.'_current"> Only From Current Site</label></div>';
					$echo_str .= $dlp."\t\t".'</td></tr><tr>';
				}
				$echo_str .= $dlp."\t\t".'<td><select id="wfu_attribute_'.$attr.'" multiple="multiple" size="'.min(count($blogs), 10).'" onchange="wfu_update_bloglist_value(\''.$attr.'\');"'.( $only_current || in_array('all', $selected) ? ' disabled="disabled"' : '' ).'>';
				foreach ( $blogs as $blog )
					$echo_str .= $dlp."\t\t\t".'<option value="'.$blog->blog_id.'"'.( in_array($blog->blog_id, $blogsselected) ? ' selected="selected"' : '' ).'>'.$blog->path.'</option>';
				$echo_str .= $dlp."\t\t".'</select>';
				$echo_str .= $dlp."\t\t".'</td><td>';
				$echo_str .= $dlp."\t\t".'<div class="wfu_bloglist_checkbtn"><input id="wfu_attribute_'.$attr.'_all" type="checkbox" onchange="wfu_update_bloglist_value(\''.$attr.'\');"'.( in_array('all', $selected) ? ' checked="checked"' : '' ).( $only_current ? ' disabled="disabled"' : '' ).' /><label for="wfu_attribute_'.$attr.'_all"> Select all</label></div>';
				$echo_str .= $dlp."\t\t".'</td></tr></tbody></table>';
			}
			else {
				$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
				$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="text" name="wfu_text_elements" value="'.$val.'" />';				
			}
		}
		elseif ( $def['type'] == "stringmatch" ) {
			$matchfield = "";
			$matchcriterion = "equal to";
			$matchvalue = "";
			preg_match('/^field:(.*?);\s*criterion:(.*?)\s*;\s*value:(.*)$/', $def['value'], $matches);
			if ( count($matches) == 4 ) {
				$matchfield = $matches[1];
				$matchcriterion = $matches[2];
				$matchvalue = $matches[3];
			}
			$echo_str .= $dlp."\t\t".'<div style="white-space:nowrap;">';
			$echo_str .= $dlp."\t\t".'<label>Field </label><input id="wfu_attribute_'.$attr.'_matchfield" type="text" name="wfu_stringmatch_elements" value="'.$matchfield.'" />';
			$echo_str .= $dlp."\t\t".'<select id="wfu_attribute_'.$attr.'_matchcriterion" value="'.$matchcriterion.'" onchange="wfu_update_stringmatch_value(\''.$attr.'\');">';
			$echo_str .= $dlp."\t\t\t".'<option value="equal to"'.( $matchcriterion == "equal to" ? 'selected="selected"' : '' ).'>equal to</option>';
			$echo_str .= $dlp."\t\t\t".'<option value="starts with"'.( $matchcriterion == "starts with" ? 'selected="selected"' : '' ).'>starts with</option>';
			$echo_str .= $dlp."\t\t\t".'<option value="ends with"'.( $matchcriterion == "ends with" ? 'selected="selected"' : '' ).'>ends with</option>';
			$echo_str .= $dlp."\t\t\t".'<option value="contains"'.( $matchcriterion == "contains" ? 'selected="selected"' : '' ).'>contains</option>';
			$echo_str .= $dlp."\t\t\t".'<option value="not equal to"'.( $matchcriterion == "not equal to" ? 'selected="selected"' : '' ).'>not equal to</option>';
			$echo_str .= $dlp."\t\t\t".'<option value="does not start with"'.( $matchcriterion == "does not start with" ? 'selected="selected"' : '' ).'>does not start with</option>';
			$echo_str .= $dlp."\t\t\t".'<option value="does not end with"'.( $matchcriterion == "does not end with" ? 'selected="selected"' : '' ).'>does not end with</option>';
			$echo_str .= $dlp."\t\t\t".'<option value="does not contain"'.( $matchcriterion == "does not contain" ? 'selected="selected"' : '' ).'>does not contain</option>';
			$echo_str .= $dlp."\t\t".'</select>';
			$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'_matchvalue" type="text" name="wfu_stringmatch_elements" value="'.$matchvalue.'" />';
			$echo_str .= $dlp."\t\t".'</div>';
		}
		elseif ( $def['type'] == "columns" ) {
			$selected = explode(",", $def['value']);
			if ( count($selected) == 1 && $selected[0] == "" ) $selected = array();
			$echo_str .= $dlp."\t".'<table class="wfu_columns_container"><tbody><tr>';
			$echo_str .= $dlp."\t\t".'<td>';
			$echo_str .= $dlp."\t\t\t".'<label class="wfu_columns_listtitle">Available Columns</label><br />';
			$echo_str .= $dlp."\t\t\t".'<select id="wfu_attribute_'.$attr.'_sourcelist" multiple="multiple" size="'.min(count($def['listitems']), 10).'">';
			foreach ( $def['listitems'] as $item ) {
				$item_required = ( substr($item, 0, 1) == "*" );
				if ( $item_required ) $item = substr($item, 1);
				$item_parts = explode("/", $item, 2);
				if ( count($item_parts) == 1 ) {
					$item_name = $item;
					$item_label = $item.( $item_required ? ' (required)' : '' );
				}
				else {
					$item_name = $item_parts[0];
					$item_label = $item_parts[1];
				}
				$itemprops[$item_name] = array( 'label' => $item_label, 'required' => $item_required );
				$echo_str .= $dlp."\t\t\t\t".'<option value="'.$item_name.'"'.( $item_required ? ' class="wfu_columns_item_required"' : '' ).' onclick="wfu_columns_itemclicked(this);">'.$item_label.'</option>';
			}
			$echo_str .= $dlp."\t\t\t".'</select>';
			$echo_str .= $dlp."\t\t".'</td>';
			foreach ( $itemprops as $item_name => $prop )
				if ( $prop['required'] && !in_array($item_name, $selected) )
					array_splice($selected, 0, 0, array( $item_name ));
			$echo_str .= $dlp."\t\t".'<td align="center">';
			$echo_str .= $dlp."\t\t\t".'<label class="wfu_columns_listtitle"></label><br />';
			$echo_str .= $dlp."\t\t\t".'<button class="wfu_columns_addbutton" onclick="wfu_columns_buttonaction(\''.$attr.'\', \'add\');">Add &gt;&gt;</button><br />';
			$echo_str .= $dlp."\t\t\t".'<button class="wfu_columns_addbutton" onclick="wfu_columns_buttonaction(\''.$attr.'\', \'del\');">&lt;&lt; Remove</button>';
			$echo_str .= $dlp."\t\t".'</td>';
			$echo_str .= $dlp."\t\t".'<td>';
			$echo_str .= $dlp."\t\t\t".'<label class="wfu_columns_listtitle">Displayed Columns</label><br />';
			$echo_str .= $dlp."\t\t\t".'<select id="wfu_attribute_'.$attr.'" multiple="multiple" size="'.min(count($def['listitems']), 10).'" onchange="wfu_update_columns(\''.$attr.'\');">';
			foreach ( $selected as $item_name ) {
				$echo_str .= $dlp."\t\t\t\t".'<option value="'.$item_name.'"'.( $itemprops[$item_name]['required'] ? ' class="wfu_columns_item_required"' : '' ).' onclick="wfu_columns_itemclicked(this);">'.$itemprops[$item_name]['label'].'</option>';
			}
			$echo_str .= $dlp."\t\t\t".'</select>';
			$echo_str .= $dlp."\t\t".'</td>';
			$echo_str .= $dlp."\t\t".'<td align="center">';
			$echo_str .= $dlp."\t\t\t".'<label class="wfu_columns_listtitle"></label><br />';
			$echo_str .= $dlp."\t\t\t".'<button class="wfu_columns_addbutton" onclick="wfu_columns_buttonaction(\''.$attr.'\', \'up\');">&#8593;</button><br />';
			$echo_str .= $dlp."\t\t\t".'<button class="wfu_columns_addbutton" onclick="wfu_columns_buttonaction(\''.$attr.'\', \'down\');">&#8595;</button>';
			$echo_str .= $dlp."\t\t".'</td>';
			$echo_str .= $dlp."\t".'</tr></tbody></table>';
			
		}
		elseif ( $def['type'] == "dimensions" ) {
			$vals_arr = explode(",", $def['value']);
			$vals = array();
			foreach ( $vals_arr as $val_raw ) {
				if ( trim($val_raw) != "" ) {
					list($val_id, $val) = explode(":", $val_raw);
					$vals[trim($val_id)] = trim($val);
				}
			}
			$dims = array();
			foreach ( $components as $comp ) {
				if ( $comp['dimensions'] == null ) $dims[$comp['id']] = $comp['name'];
				else foreach ( $comp['dimensions'] as $dimraw ) {
					list($dim_id, $dim_name) = explode("/", $dimraw);
					$dims[$dim_id] = $dim_name;
				}
			}
			foreach ( $dims as $dim_id => $dim_name ) {
				if ( !array_key_exists($dim_id, $vals) ) $vals[$dim_id] = "";
				$echo_str .= $dlp."\t\t".'<span style="display:inline-block; width:130px;">'.$dim_name.'</span><input id="wfu_attribute_'.$attr.'_'.$dim_id.'" type="text" name="wfu_dimension_elements_'.$attr.'" class="wfu_short_text" value="'.$vals[$dim_id].'" /><br />';
			}
		}
		elseif ( $def['type'] == "userfields" ) {
			$fields_arr = explode("/", $def['value']);
			$fields = array();
			foreach ( $fields_arr as $field_raw ) {
				$is_req = ( substr($field_raw, 0, 1) == "*" );
				if ( $is_req ) $field_raw = substr($field_raw, 1);
				if ( $field_raw != "" ) array_push($fields, array( "name" => $field_raw, "required" => $is_req ));
			}
			if ( count($fields) == 0 ) array_push($fields, array( "name" => "", "required" => false ));
			$echo_str .= $dlp."\t\t".'<div id="wfu_attribute_'.$attr.'" class="wfu_userdata_container">';
			foreach ( $fields as $field ) {
				$echo_str .= $dlp."\t\t\t".'<div class="wfu_userdata_line">';
				$echo_str .= $dlp."\t\t\t\t".'<input type="text" name="wfu_userfield_elements" value="'.$field['name'].'" />';
				$echo_str .= $dlp."\t\t\t\t".'<div class="wfu_userdata_action" onclick="wfu_userdata_add_field(this);"><img src="'.WFU_IMAGE_ADMIN_USERDATA_ADD.'" ></div>';
				$echo_str .= $dlp."\t\t\t\t".'<div class="wfu_userdata_action wfu_userdata_action_disabled" onclick="wfu_userdata_remove_field(this);"><img src="'.WFU_IMAGE_ADMIN_USERDATA_REMOVE.'" ></div>';
				$echo_str .= $dlp."\t\t\t\t".'<input type="checkbox"'.( $field['required'] ? 'checked="checked"' : '' ).' onchange="wfu_update_userfield_value({target:this});" />';
				$echo_str .= $dlp."\t\t\t\t".'<span>Required</span>';
				$echo_str .= $dlp."\t\t\t".'</div>';
			}
			$echo_str .= $dlp."\t\t".'</div>';
		}
		elseif ( $def['type'] == "color" ) {
			$val = str_replace(array( "%n%", "%dq%", "%brl%", "%brr%" ), array( "\n", "&quot;", "[", "]" ), $def['value']);
			$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="text" name="wfu_text_elements" class="wfu_color_field" value="'.$val.'" />';
		}
		elseif ( $def['type'] == "color-triplet" ) {
			$triplet = explode(",", $def['value']);
			foreach ( $triplet as $key => $item ) $triplet[$key] = trim($item);
			if ( count($triplet) == 2 ) $triplet = array( $triplet[0], $triplet[1], "#000000");
			elseif ( count($triplet) == 1 ) $triplet = array( $triplet[0], "#FFFFFF", "#000000");
			elseif ( count($triplet) < 3 ) $triplet = array( "#000000", "#FFFFFF", "#000000");
			$echo_str .= $dlp."\t\t".'<div class="wfu_color_container"><label style="display:inline-block; width:120px; margin-top:-16px;">Text Color</label><input id="wfu_attribute_'.$attr.'_color" type="text" class="wfu_color_field" name="wfu_triplecolor_elements" value="'.$triplet[0].'" /></div>';
			$echo_str .= $dlp."\t\t".'<div class="wfu_color_container"><label style="display:inline-block; width:120px; margin-top:-16px;">Background Color</label><input id="wfu_attribute_'.$attr.'_bgcolor" type="text" class="wfu_color_field" name="wfu_triplecolor_elements" value="'.$triplet[1].'" /></div>';
			$echo_str .= $dlp."\t\t".'<div class="wfu_color_container"><label style="display:inline-block; width:120px; margin-top:-16px;">Border Color</label><input id="wfu_attribute_'.$attr.'_borcolor" type="text" class="wfu_color_field" name="wfu_triplecolor_elements" value="'.$triplet[2].'" /></div>';
		}
		else {
			$echo_str .= $dlp."\t\t".'<input id="wfu_attribute_'.$attr.'" type="text" name="wfu_text_elements" value="'.$def['value'].'" />';
			if ( $def['variables'] != null ) $echo_str .= $dlp.wfu_insert_variables($def['variables'], 'wfu_variable wfu_variable_'.$attr);
		}
		$echo_str .= $dlp."\t".'</div></td>';
		if ( $def['parent'] == "" ) {
			$echo_str .= $dlp."\t".'<td style="position:relative; vertical-align:top; padding:0;"><div class="wfu_td_div">';
			$block_open = false;
		}
		else {
			$echo_str .= $dlp.'</tr>';
			$subblock_open = true;						
		}
	}
	if ( $subblock_open ) {
		$echo_str .= "\n\t\t\t\t\t\t".'</div>';
	}
	if ( $block_open ) {
		$echo_str .= "\n\t\t\t\t\t".'</div></td>';
		$echo_str .= "\n\t\t\t\t".'</tr>';
	}
	$echo_str .= "\n\t\t\t".'</tbody>';
	$echo_str .= "\n\t\t".'</table>';
	$echo_str .= "\n\t".'</div>';
	$echo_str .= "\n\t".'<div id="wfu_global_dialog_container" class="wfu_global_dialog_container">';
	$echo_str .= "\n\t".'</div>';
	$handler = 'function() { wfu_Attach_Admin_Events('.( $data == "" ? 'true' : 'false' ).'); }';
	$echo_str .= "\n\t".'<script type="text/javascript">if(window.addEventListener) { window.addEventListener("load", '.$handler.', false); } else if(window.attachEvent) { window.attachEvent("onload", '.$handler.'); } else { window["onload"] = '.$handler.'; }</script>';
	$echo_str .= "\n".'</div>';
//	$echo_str .= "\n\t".'<div style="margin-top:10px;">';
//	$echo_str .= "\n\t\t".'<label>Final shortcode text</label>';
//	$echo_str .= "\n\t".'</div>';

	echo $echo_str;
}

function wfu_insert_variables($variables, $class) {
	$ret = "";
	foreach ( $variables as $variable )
		if ( $variable == "%userdataXXX%" ) $ret .= "\t\t".'<select class="'.$class.'" name="wfu_userfield_select" title="'.constant("WFU_VARIABLE_TITLE_".strtoupper(str_replace("%", "", $variable))).'" onchange="wfu_insert_userfield_variable(this);"><option style="display:none;">%userdataXXX%</option></select>';
		elseif ( $variable != "%n%" && $variable != "%dq%" && $variable != "%brl%" && $variable != "%brr%" ) $ret .= "\t\t".'<span class="'.$class.'" title="'.constant("WFU_VARIABLE_TITLE_".strtoupper(str_replace("%", "", $variable))).'" ondblclick="wfu_insert_variable(this);">'.$variable.'</span>';
	return $ret;
}

?>
